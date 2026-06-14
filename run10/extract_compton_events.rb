#!/usr/bin/env ruby

require 'comptonsoft'
require 'fileutils'
require 'csv'

class MyApp < ANL::ANLApp
  attr_accessor :inputs, :output

  def setup()
    add_namespace ComptonSoft

    chain :CSHitCollection
    chain :ConstructDetector
    with_parameters(detector_configuration: "database/detector_configuration.xml",
                    verbose_level: 1)
    chain :ReadHitTree
    with_parameters(file_list: @inputs)
    chain :EventReconstruction
    with_parameters(reconstruction_method: "NanoGRAMS",
                    source_distant: false,
                    source_position: vec(50.0, 0.0, 0.0),
                    parameter_file: "parfile_NanoGRAMS.yaml")
    chain :WriteComptonEventTree
    chain :SaveData
    with_parameters(output: @output)
  end
end


### main ###
data_type = "HSTD14"
filename      = "metadata/run10data_with_interpolated_FEC_#{data_type}.csv"
outdir_parent = "products"
outdir = outdir_parent
hittree_files = []

CSV.foreach(filename, headers: true) do |row|
  time_array = row["time"].split("/")
  hittree_path = "#{outdir_parent}/#{time_array[0]}/#{time_array[1]}/hittree.root"
  hittree_files << hittree_path
end

### Event reconstruction
a = MyApp.new
a.console = false
a.inputs = hittree_files
a.output = "#{outdir}/compton_#{data_type}.root"

a.run(:all)
