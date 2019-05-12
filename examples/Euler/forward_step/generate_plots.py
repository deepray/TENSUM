import sys

if((len(sys.argv) > 3 and len(sys.argv) < 6) or len(sys.argv) < 3):
    print('ERROR: Must use the script as follows')
    print('   visit -cli -s generate_plots <SAMPLE_DIR> <variable>')
    print('or')
    print('   visit -cli -s generate_plots <SAMPLE_DIR> <variable> <cmin> <cmax> <N>')
    print('where:')
    print('   <SAMPLE_DIR>  --> Directory containing solution files')
    print('   <variable>    --> a valid scalar variable name')
    print('   <cmin>        --> minimum contour and pseudocolor level')
    print('   <cmax>        --> maximum contour and pseudocolor level')
    print('   <N>           --> number contour lines')
    exit()

DB_Dir      = sys.argv[1]
var         = sys.argv[2]

    
master_file = DB_Dir+'/master_file.visit'

OpenDatabase(master_file)

DefineScalarExpression("velx","velocity[0]")
DefineScalarExpression("vely","velocity[1]")

# Draw plots
AddPlot("Pseudocolor", var)
AddPlot("Contour",var)
DrawPlots()

# PseudoColor attributes
PCAtts = PseudocolorAttributes()
PCAtts.colorTableName = "hot_desaturated"
if(len(sys.argv) > 3):
    PCAtts.minFlag = 1
    PCAtts.maxFlag = 1
    PCAtts.min = float(sys.argv[3])
    PCAtts.max = float(sys.argv[4])
SetPlotOptions(PCAtts)

# Contour attributes
ContourAtts = ContourAttributes()
ContourAtts.colorType=ContourAtts.ColorBySingleColor
ContourAtts.singleColor=(0,0,0,255)
ContourAtts.legendFlag = 0
ContourAtts.lineWidth = 1
if(len(sys.argv) > 3):
    ContourAtts.minFlag = 1
    ContourAtts.maxFlag = 1
    ContourAtts.min = float(sys.argv[3])
    ContourAtts.max = float(sys.argv[4])
    ContourAtts.contourNLevels = int(sys.argv[5])
else:
    ContourAtts.contourNLevels=20
SetPlotOptions(ContourAtts)

Nt = TimeSliderGetNStates()

# Save attributes
s = SaveWindowAttributes()
s.family = 0
s.outputToCurrentDirectory = 0
s.outputDirectory = "."
SetSaveWindowAttributes(s)


# Changing annotations
p = AnnotationAttributes()
p.axes2D.visible = 0
p.userInfoFlag = 0
p.databaseInfoFlag = 0
SetAnnotationAttributes(p)

# Zoom into airfoil
MoveAndResizeWindow(1,0,0,700,300)
# v = GetView2D()
# #v.windowCoords = (-3, 3,-3,3)
# v.viewportCoords = (-0.5,1,0,0.7)
# SetView2D(v)

# Editting legend
plotName = GetPlotList().GetPlots(0).plotName
legend = GetAnnotationObject(plotName)
legend.drawMinMax = 1
legend.drawTitle = 0
legend.managePosition = 0
legend.xScale = 1
legend.yScale = 3
legend.fontHeight = 0.05
legend.fontBold = 1
legend.position = (0.01,0.85)
legend.numberFormat = "%# -6.2e"

for t in range(Nt):

    TimeSliderSetState(t)
    s.format=s.PNG
    s.fileName="%s_%d.png"%(var,t)
    SetSaveWindowAttributes(s)
    SaveWindow()

exit()
