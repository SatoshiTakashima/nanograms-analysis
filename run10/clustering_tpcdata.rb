#!/usr/bin/env ruby

require 'comptonsoft'
require 'fileutils'
require 'csv'

class MyAppDataReduction < ANL::ANLApp
  attr_accessor :tpc_tree_file, :hittree_file
  attr_accessor :rawhittree_file
  attr_accessor :quicklook_file
  attr_accessor :gain_tp_file, :gain_tp_hash

  def setup
    add_namespace ComptonSoft

    chain :NanoGRAMSHitExtraction
    with_parameters(
      config_file:     "config_pipeline.yaml",
      tpctree_file:    @tpc_tree_file,
      rawhittree_file: @rawhittree_file,
      quicklook_file:  @quicklook_file,
    )

    chain :NanoGRAMSCalibration
    calibration_parameters = { hittree_file: @hittree_file }
    calibration_parameters[:gain_tp_hash] = @gain_tp_hash if @gain_tp_hash
    calibration_parameters[:gain_tp_file] = @gain_tp_file if @gain_tp_file
    with_parameters(**calibration_parameters)
  end
end

def merge_root_files(output_file, input_files)
  existing_files = input_files.select { |path| File.file?(path) }
  missing_files = input_files - existing_files

  missing_files.each do |path|
    warn("[merge] skip missing file: #{path}")
  end

  if existing_files.empty?
    warn("[merge] no input files for #{output_file}")
    return
  end

  FileUtils.mkdir_p(File.dirname(output_file))
  puts("[merge] #{output_file}")
  unless system("hadd", "-f", output_file, *existing_files)
    raise "hadd failed: #{output_file}"
  end
end


### main ###
filename      = "metadata/run10data_with_interpolated_FEC_HSTD14.csv"
#filename      = "metadata/run10data_with_interpolated_FEC.csv"
data_root_dir = "/Users/takashima/work/grams/run/run10/data/tpc/data"
gain_tp_file  = "testpulse_analysis/products/run10_testpulse_data.csv"
outdir_parent = "products"
hittree_files = []
compton_files = []
CSV.foreach(filename, headers: true) do |row|
  time_array = row["time"].split("/")

 data_dir = "#{data_root_dir}/#{time_array[0]}/#{time_array[1]}"
 outdir   = "#{outdir_parent}/#{time_array[0]}/#{time_array[1]}"
 FileUtils.mkdir_p(outdir)
 
 ### Data reduction
 a = MyAppDataReduction.new
 a.console = false
 a.tpc_tree_file   = "#{data_dir}/tpc_data.root"
 a.rawhittree_file = "#{outdir}/rawhittree.root"
 a.quicklook_file  = "#{outdir}/quicklook_tree.root"
 a.gain_tp_file    = gain_tp_file
 a.hittree_file    = "#{outdir}/hittree.root"


 a.run(:all)
 hittree_files << a.hittree_file
end

merge_root_files("#{outdir_parent}/hittree_merged.root", hittree_files)
