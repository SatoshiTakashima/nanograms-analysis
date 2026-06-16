#!/usr/bin/env ruby

require 'comptonsoft'
require 'fileutils'
require 'csv'

class MyApp < ANL::ANLApp
  attr_accessor :inputs, :output
  attr_accessor :source_position

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
                    source_position: @source_position,
                    parameter_file: "parfile_NanoGRAMS.yaml")
    chain :WriteComptonEventTree
    chain :SaveData
    with_parameters(output: @output)
  end
end


### main ###
#data_type     = "HSTD14"
data_type     = "zm4cm"
source_position = vec(50.0, 0.0, -4.0)
filename      = "metadata/data_group_#{data_type}.csv"
outdir_parent = "products"
outdir        = "#{outdir_parent}/#{data_type}"
hittree_files = []

FileUtils.mkdir_p(outdir)

CSV.foreach(filename, headers: true) do |row|
  time_array = row["time"].split("/")
  hittree_path = "#{outdir_parent}/#{time_array[0]}/#{time_array[1]}/hittree.root"
  hittree_files << hittree_path
end

### Event reconstruction
a = MyApp.new
a.console = false
a.source_position = source_position
a.inputs = hittree_files
a.output = "#{outdir}/compton.root"

a.run(:all)
