#!/usr/bin/env ruby

require 'comptonsoft'
require 'fileutils'
require 'csv'
require 'open3'

Dir.chdir(File.expand_path(__dir__))

class MyAppDataReduction < ANL::ANLApp
  attr_accessor :tpc_tree_file, :hittree_file
  attr_accessor :quicklook_file
  attr_accessor :gain_tp_file, :gain_tp_hash
  attr_accessor :gain_cache_seconds

  def setup
    add_namespace ComptonSoft

    chain :NanoGRAMSHitExtraction
    extraction_parameters = {
      config_file:     "metadata/config_pipeline.yaml",
      tpctree_file:    @tpc_tree_file,
    }
    extraction_parameters[:quicklook_file] = @quicklook_file if @quicklook_file
    with_parameters(**extraction_parameters)

    chain :NanoGRAMSCalibration
    calibration_parameters = { hittree_file: @hittree_file }
    calibration_parameters[:gain_tp_hash] = @gain_tp_hash if @gain_tp_hash
    calibration_parameters[:gain_tp_file] = @gain_tp_file if @gain_tp_file
    calibration_parameters[:gain_cache_seconds] = @gain_cache_seconds if @gain_cache_seconds
    with_parameters(**calibration_parameters)
  end
end

class MyAppMergeHitTree < ANL::ANLApp
  attr_accessor :hittree_files, :output

  def setup
    add_namespace ComptonSoft

    chain :MergeHitTree
    with_parameters(
      hittree_files: @hittree_files,
      output:        @output
    )
  end
end


### main ###
data_root_dir = "/Users/grams/run/run12/tpc_data"
outdir_parent = "products/eventfile"
gain_cache_seconds = 60.0
python = ENV.fetch("PYTHON", "/Users/grams/.pyenv/versions/3.14.6/bin/python")
quicklook_enabled = !%w[0 false no off].include?(ENV.fetch("MAKE_QUICKLOOK", "1").downcase)
quicklook_chunk_entries = Integer(ENV.fetch("QUICKLOOK_CHUNK_ENTRIES", "4000"))

def event_dir_name(time_id)
  time_id.strip.tr("/", "_")
end

def tpc_data_dir(data_root_dir, time_id)
  flat_dir = File.join(data_root_dir, event_dir_name(time_id))
  return flat_dir if Dir.exist?(flat_dir)

  File.join(data_root_dir, *time_id.strip.split("/"))
end

def root_file_has_tree?(file_path, tree_name)
  output, status = Open3.capture2("rootls", file_path)
  status.success? && output.lines.map(&:strip).include?(tree_name)
end

def comptonsoft_time_id(time_id)
  text = time_id.strip
  match = text.match(/\A(\d{8})_(\d{4}_\d{2})\z/)
  return "#{match[1]}/#{match[2]}" if match

  text
end

def write_comptonsoft_gain_csv(source, output)
  table = CSV.read(source, headers: true)
  CSV.open(output, "w") do |csv|
    csv << table.headers
    table.each do |row|
      row["time_id"] = comptonsoft_time_id(row["time_id"]) if row["time_id"]
      csv << row
    end
  end
end

def split_tpctree_for_quicklook(python, file_path, chunk_dir, max_entries)
  config_path = File.join(File.dirname(file_path), "config_dpp.yaml")
  command = [
    python,
    File.expand_path("tools/split_tpctree.py", __dir__),
    file_path,
    chunk_dir,
    "--max-entries", max_entries.to_s,
    "--config", config_path,
  ]
  output, status = Open3.capture2e(*command)
  raise "Failed to split #{file_path}:\n#{output}" unless status.success?

  output.lines.filter_map do |line|
    part, path = line.strip.split("\t", 2)
    [part, path] if part && path
  end
end

data_group_list = ["test"]
only_time_id = ENV["TIME_ID"]&.strip

data_group_list.each do |tag|
  source_gain_tp_file = File.expand_path("products/testpulse_data.csv", __dir__)
  gain_tp_file  = File.expand_path("products/testpulse_data_comptonsoft.csv", __dir__)
  write_comptonsoft_gain_csv(source_gain_tp_file, gain_tp_file)
  filename      = "metadata/data_group/data_group_#{tag}.csv"
  hittree_files = []
  compton_files = []
  CSV.foreach(filename, headers: true) do |row|
    time_id = row["time"]
    next if only_time_id && event_dir_name(time_id) != event_dir_name(only_time_id)

    data_dir = tpc_data_dir(data_root_dir, time_id)
    outdir   = File.join(outdir_parent, event_dir_name(time_id))
    FileUtils.mkdir_p(outdir)

    hittree_local_files = []
    filePattern = "tpc_data*.root"
    Dir.glob(File.join(data_dir, filePattern)).sort.each do |file_path|
      unless root_file_has_tree?(file_path, "tpctree")
        warn "[skip] missing tpctree: #{file_path}"
        next
      end

      basename = File.basename(file_path)
      match = basename.match(/tpc_data_?(\d+)\.root/)
      next unless match
      num_str = match[1]

      inputs = if quicklook_enabled
                 chunk_dir = File.join(outdir, "tpctree_chunks", "tpc_data#{num_str}")
                 split_tpctree_for_quicklook(python, file_path, chunk_dir, quicklook_chunk_entries)
               else
                 [[nil, file_path]]
               end

      inputs.each do |part, input_file_path|
        suffix = part ? "#{num_str}_part#{part}" : num_str
        hittree_file = "#{outdir}/hittree#{suffix}.root"
        quicklooktree_file = quicklook_enabled ? "#{outdir}/quicklook#{suffix}.root" : nil
        a = MyAppDataReduction.new
        a.console = false
        a.tpc_tree_file = input_file_path
        a.gain_tp_file = gain_tp_file
        a.gain_cache_seconds = gain_cache_seconds
        a.hittree_file = hittree_file
        a.quicklook_file = quicklooktree_file if quicklooktree_file
        a.run(:all)
        hittree_local_files << hittree_file
      end
    end

    #### Merging hittree
    next if hittree_local_files.empty?
    b = MyAppMergeHitTree.new
    b.console = false
    b.hittree_files = hittree_local_files
    b.output = "#{outdir}/hittree.root"
    b.run(:all)
    hittree_files << b.output
  end
end
