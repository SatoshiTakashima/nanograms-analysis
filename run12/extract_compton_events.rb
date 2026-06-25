#!/usr/bin/env ruby

require 'comptonsoft'
require 'fileutils'
require 'csv'

Dir.chdir(File.expand_path(__dir__))

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
                    parameter_file: "metadata/parfile_NanoGRAMS.yaml")
    chain :WriteComptonEventTree
    chain :SaveData
    with_parameters(output: @output)
  end
end


### main ###
source_position = vec(50.0, 0.0, 0.0)
outdir_parent = "products/eventfile"

def event_dir_name(time_id)
  time_id.strip.tr("/", "_")
end

data_group_list = ["test"]

data_group_list.each do |tag|
  filename = "metadata/data_group/data_group_#{tag}.csv"

  CSV.foreach(filename, headers: true) do |row|
    outdir = File.join(outdir_parent, event_dir_name(row["time"]))
    FileUtils.mkdir_p(outdir)

    ### Event reconstruction
    a = MyApp.new
    a.console = false
    a.source_position = source_position
    a.inputs = ["#{outdir}/hittree.root"]
    a.output = "#{outdir}/compton.root"

    a.run(:all)
  end
end
