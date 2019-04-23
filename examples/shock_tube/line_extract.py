# Script to extract solution along a line from partitioned vtk files

import sys
import math
import vtk
import numpy as np

# Check if a vtk file argument is given
if len(sys.argv) < 3:
   print "Please specify a vtk file name (base) followed by the number of partitioned files."
   print "Example:"
   print "  ", sys.argv[0], " soln_0003 10"
   print "will extract line data from soln_0003_0.vtk to soln_0003_10.vtk"
   sys.exit(1)

Npart = int(sys.argv[2])

# Specify line parameters: data extracted between points (x0,y) and (xN,y)
x0 = 0.0;
xN = 1.0;
y  = 0.0;
z  = 0.0;

x   = []
rho = []
ux  = []
pre = [] 

for n in range(Npart):
   # Set the vtk input file name
   vtkfile  = sys.argv[1]
   vtkfile += "_" + str(n) + ".vtk"

   # Read the unstructured grid data
   reader = vtk.vtkUnstructuredGridReader()
   reader.ReadAllScalarsOn()
   reader.ReadAllVectorsOn()
   reader.SetFileName(vtkfile)
   reader.Update()

   # Create the line source to use for the probe lines.
   line = vtk.vtkLineSource()
   line.SetPoint1(x0,y,z)
   line.SetPoint2(xN,y,z)
   line.SetResolution(1000)

   # Move the line into place and create the probe filter. For
   # vtkProbeFilter, the probe line is the input, and the underlying data
   # set is the source.
   probe = vtk.vtkProbeFilter()
   probe.SetInputConnection(line.GetOutputPort())
   probe.SetSource(reader.GetOutput())
   probe.Update()
   data=probe.GetOutput()

   # Extract velocity, density and pressure from point data
   ptdata         = data.GetPointData()
   arrayid        = ptdata.SetActiveVectors("velocity")
   velocity       = ptdata.GetArray(arrayid)
   valid_point_id = data.GetPointData().SetActiveScalars("vtkValidPointMask")
   valid_pts      = ptdata.GetArray(valid_point_id)
   arrayid        = ptdata.SetActiveScalars("density")
   density        = ptdata.GetArray(arrayid)
   arrayid        = ptdata.SetActiveScalars("pressure")
   pressure       = ptdata.GetArray(arrayid)


   for i in range(velocity.GetNumberOfTuples()):
      if valid_pts.GetTuple1(i) == 1.0:
         pt  = data.GetPoint(i)
         if x.count(pt[0]) == 0:
            x.append(pt[0])
            a  = velocity.GetTuple3(i)
            d  = density.GetTuple1(i)
            p  = pressure.GetTuple1(i)
            rho.append(d)
            ux.append(a[0])
            pre.append(p)

# Writing data to file
f = open('x.dat','w') 
idx = np.argsort(x) 

for i in range(len(idx)):
   j = idx[i]
   s  = str(x[j]) + " " + str(rho[j]) + " " + str(ux[j]) + " " + str(pre[j]) + "\n"              
   f.write(s)  

f.close()    

