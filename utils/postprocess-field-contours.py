"""
Plot over line for fields from several solutions in a L0-parametrized sequence. 
Investigate Quasistationary Transition to Kerr BH.

Usage:

Save to solution_sequence/adaptivesolver/visualization/ 
Run using paraview python shell.

"""

from paraview.simple import *
import os, glob


#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Parameters
field = 'NU'
solution_sequence = 'tss_L08'
file_name = solution_sequence + 'test' + '_%s' % (field)


# Set directories (must be hardcoded)
solution_dir = '/Users/elleryames/Research/projects/projectgeco/geon/paper_data/tss_L08/adaptive_solver'

# Get list of solutions/files
steps = [60, 66, 94] #[5, 25, 50, 60, 65, 66, 70, 80, 94]
#step_dir = glob.glob(solution_dir + '/step_')
files = [solution_dir + '/step_0{:02}/{:}_100.xdmf'.format(s, field) for s in steps]


# Set up View
# Create a new 'Line Chart View'
CreateLayout('line-source-traces')
traceView = CreateXYPlotView()
#traceView.ViewSize = [1200, 800]

def get_line_source(source):
    # create a new 'Plot Over Line'
    plotOverLine = PlotOverLine(Input=source,
    Source='High Resolution Line Source')

    # Properties modified on plotOverLine1
    plotOverLine.Tolerance = 2.22044604925031e-16

    # Properties modified on plotOverLine1.Source
    plotOverLine.Source.Point2 = [10.0, 0.0, 0.0]
    plotOverLine.Source.Resolution = 1000

    return plotOverLine

def get_data_name(xdmf_file):

    with open(xdmf_file, 'r') as infofile:
        data=infofile.read()
        items = data.split('<')
        for i in items:
            if 'Attribute Name=' in i:
                data_name = i.split('=')[1].split('"')[1]
            
    return data_name
    

# create a new 'XDMF Reader'
cval = [ [0.889998, 0.500008, 0.110002], [0.0, 0.0, 0.0], [0.0, 0.0, 0.9] ]
s = 0
CreateLayout('renders')
for f, step, c in zip(files, steps, cval):
    #load data
    reader = XDMFReader(FileNames=f)
    data_name = get_data_name(f)

    # render solution
    #CreateLayout('step %s' % s)
    renderView = CreateRenderView('RenderView')
    #renderView.ViewSize = [600, 400]
    display = Show(reader, renderView)
    s += 1

    # take trace and plot
    trace = get_line_source(reader)
    plotOverLineDisplay = Show(trace, traceView)
    # Edit legend name
    plotOverLineDisplay.SeriesLabel = ['arc_length', 'arc_length', data_name, 'step {:}'.format(step), 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
    # Edit line color
    plotOverLineDisplay.SeriesColor = ['arc_length', '0', '0', '0', data_name, '{:}'.format(c[0]), '{:}'.format(c[1]), '{:}'.format(c[2]), 'vtkValidPointMask', '0.220005', '0.489998', '0.719997', 'Points_X', '0.300008', '0.689998', '0.289998', 'Points_Y', '0.6', '0.310002', '0.639994', 'Points_Z', '1', '0.500008', '0', 'Points_Magnitude', '0.650004', '0.340002', '0.160006']
    Hide(reader, renderView)

# Save screenshot
GetLayoutByName('line-source-traces')
SaveScreenshot(os.path.join(solution_dir + '/visualization/contour_plots/' + file_name + ".png"), traceView)  

   
