#!/usr/bin/env ruby

require 'comptonsoft'
require 'fileutils'
require 'csv'

Dir.chdir(File.expand_path(__dir__))

class MyAppDataReduction < ANL::ANLApp
  attr_accessor :tpc_tree_file, :hittree_file
  attr_accessor :quicklook_file
  attr_accessor :gain_tp_file, :gain_tp_hash

  def setup
    add_namespace ComptonSoft

    chain :NanoGRAMSHitExtraction
    with_parameters(
      config_file:     "metadata/config_pipeline.yaml",
      tpctree_file:    @tpc_tree_file,
      quicklook_file:  @quicklook_file,
    )

    chain :NanoGRAMSCalibration
    calibration_parameters = { hittree_file: @hittree_file }
    calibration_parameters[:gain_tp_hash] = @gain_tp_hash if @gain_tp_hash
    calibration_parameters[:gain_tp_file] = @gain_tp_file if @gain_tp_file
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
data_root_dir = "/Users/takashima/work/grams/run/run12/mockdata/tpc/data"
outdir_parent = "products"

data_group_list = ["test"]

data_group_list.each do |tag|
  gain_tp_file  = "../products/interpolated_gain_#{tag}.csv"
  filename      = "metadata/data_group/data_group_#{tag}.csv"
  hittree_files = []
  compton_files = []
  CSV.foreach(filename, headers: true) do |row|
    time_array = row["time"].split("/")
  
    data_dir = "#{data_root_dir}/#{time_array[0]}/#{time_array[1]}"
    outdir   = "#{outdir_parent}/#{time_array[0]}/#{time_array[1]}"
    FileUtils.mkdir_p(outdir)

    hittree_local_files = []
    filePattern = "tpc_data[0-9][0-9][0-9][0-9][0-9].root"
    Dir.glob(File.join(data_dir, filePattern)).sort.each do |file_path|

      basename = File.basename(file_path)
      match = basename.match(/tpc_data(\d{5})\.root/)
      next unless match
      num_str = match[1]
      p file_path

      quicklook_file = "#{outdir}/quicklook_tree#{num_str}.root"
      hittree_file   = "#{outdir}/hittree#{num_str}.root"
      p hittree_file
      a = MyAppDataReduction.new
      a.console = false
      a.tpc_tree_file  = file_path 
      a.quicklook_file = quicklook_file
      a.gain_tp_file   = gain_tp_file
      a.hittree_file   = hittree_file
      a.run(:all)
      hittree_local_files << hittree_file
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
