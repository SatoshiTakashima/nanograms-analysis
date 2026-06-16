#!/usr/bin/env ruby

require 'comptonsoft'

Degree = Math::PI / 180.0

class MyApp < ANL::ANLApp
  attr_accessor :inputs, :output

  def setup()
    add_namespace ComptonSoft

    chain :ConstructDetector
    with_parameters(detector_configuration: "database/detector_configuration.xml")
    chain :ReadComptonEventTree
    with_parameters(file_list: inputs)
    chain :BackProjection
    with_parameters(x_min: -20.0,
                    x_max: +20.0,
                    y_min: -20.0,
                    y_max: +20.0,
                    plane_normal: vec(1.0, 0.0, 0.0),
                    plane_point: vec(50.0, 0.0, 0.0),
                    plane_yaxis:  vec(0.0, 0.0, 1.0),
                    num_pixel_x: 100,
                    num_pixel_y: 100,
                    arm: 3.0,
                    num_points: 10000)
    chain :SaveData
    with_parameters(output: output)
  end
end

app = MyApp.new
app.inputs = ["products/HSTD14/compton.root"]
app.output = "products/HSTD14/compton_image.root"
app.run(:all)
