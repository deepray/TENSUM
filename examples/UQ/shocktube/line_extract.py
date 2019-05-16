import sys

Sample_start = int(sys.argv[1])
Sample_end = int(sys.argv[2])


for s in range(Sample_start,Sample_end+1):
	master_file = 'SAMPLE_'+str(s)+'/master_file.visit'
	ofile_base  = 'SAMPLE_'+str(s)+'/line_data'

	OpenDatabase(master_file)

	DefineScalarExpression("velx","velocity[0]")

	LineoutAtts = LineoutAttributes()
	LineoutAtts.point1 = (0, 0.02, 0)
	LineoutAtts.point2 = (1, 0.02, 0)
	LineoutAtts.interactive = 0
	LineoutAtts.ignoreGlobal = 0
	LineoutAtts.samplingOn = 1
	LineoutAtts.numberOfSamplePoints = 100
	LineoutAtts.reflineLabels = 0

	Nt = TimeSliderGetNStates()

	for t in range(Nt):

		TimeSliderSetState(t)

		# Density
		AddPlot("Curve", "operators/Lineout/density", 1, 1)
		SetOperatorOptions(LineoutAtts, 1)
		DrawPlots()
		linedata = GetPlotInformation()["Curve"]
		xc  = linedata[0::2]
		rho = linedata[1::2]
		ClearWindow()

		# Velocity
		AddPlot("Curve", "operators/Lineout/velx", 1, 1)
		SetOperatorOptions(LineoutAtts, 1)
		DrawPlots()
		linedata = GetPlotInformation()["Curve"]
		velx = linedata[1::2]
		ClearWindow()

		# Pressure
		AddPlot("Curve", "operators/Lineout/pressure", 1, 1)
		SetOperatorOptions(LineoutAtts, 1)
		DrawPlots()
		linedata = GetPlotInformation()["Curve"]
		pre = linedata[1::2]
		ClearWindow()

		# Write it as "x y z val"
		ofile = ofile_base + '_'+str(t)+'.dat'
		f=open(ofile,'w')
		for i in range(len(xc)): 
		   f.write("%g %g %g %g \n" % (xc[i], rho[i], velx[i], pre[i]))
		#print "%g %g" % (xc[idx], vals[idx])
		f.close()
	
	DeleteAllPlots()
	CloseDatabase(master_file)	

exit()
