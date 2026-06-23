#!/usr/bin/env ruby

require 'comptonsoft'
require 'fileutils'
require 'csv'

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


### main ###
data_root_dir = "/Users/takashima/work/grams/run/run10/data/tpc/data"
outdir_parent = "products"

data_group_list = ["z0cm", "z4cm", "zm4cm"]

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
   
   ### Data reduction
   a = MyAppDataReduction.new
   a.console = false
   a.tpc_tree_file  = "#{data_dir}/tpc_data.root"
   a.quicklook_file = "#{outdir}/quicklook_tree.root"
   a.gain_tp_file   = gain_tp_file
   a.hittree_file   = "#{outdir}/hittree.root"
  
   a.run(:all)
   hittree_files << a.hittree_file
  end
end
