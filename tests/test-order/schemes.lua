{																															-- advect wave L1-error				Sod L1-error
	{solver='weno', wenoMethod='1996 Jiang Shu', wenoFlux='Lax-Friedrichs', order=3, integrator='Runge-Kutta 4'},			-- 6.5248777651847e-08				0.0012540866319073
	{solver='weno', wenoMethod='2008 Borges', wenoFlux='Lax-Friedrichs', order=3, integrator='Runge-Kutta 4'},				-- 5.4025263260433e-09  			0.00098033101037162
	{solver='weno', wenoMethod='2010 Shen Zha', wenoFlux='Lax-Friedrichs', order=3, integrator='Runge-Kutta 4'},			-- 6.1547240085165e-07  			0.00086506093962507
	{solver='weno', wenoMethod='1996 Jiang Shu', wenoFlux='Lax-Friedrichs', order=5, integrator='Runge-Kutta 4'},			-- 2.4023578265586e-13  			0.0007741510226556
	{solver='weno', wenoMethod='2008 Borges', wenoFlux='Lax-Friedrichs', order=5, integrator='Runge-Kutta 4'},				-- 4.1535351547051e-14  			0.00064106114903531
	{solver='weno', wenoMethod='2010 Shen Zha', wenoFlux='Lax-Friedrichs', order=5, integrator='Runge-Kutta 4'},			-- 4.6202516638949e-14  			0.00061303349009959
	{solver='weno', wenoMethod='1996 Jiang Shu', wenoFlux='Roe', order=5, integrator='Runge-Kutta 4'},						-- 2.4023578265586e-13  			0.0007741510226556
	{solver='weno', wenoMethod='2008 Borges', wenoFlux='Roe', order=5, integrator='Runge-Kutta 4'},							-- 4.1535351547051e-14  			0.00064106114903531
	{solver='weno', wenoMethod='2010 Shen Zha', wenoFlux='Roe', order=5, integrator='Runge-Kutta 4'},						-- 4.6202516638949e-14  			0.00061303349009959
	{solver='weno', wenoMethod='1996 Jiang Shu', wenoFlux='Lax-Friedrichs', order=7, integrator='Runge-Kutta 4'},			-- 8.975025583835e-15   			0.00062607092123513
	{solver='weno', wenoMethod='2008 Borges', wenoFlux='Lax-Friedrichs', order=7, integrator='Runge-Kutta 4'},				-- 8.7285864300291e-15  			0.00051638820293648
	{solver='weno', wenoMethod='2010 Shen Zha', wenoFlux='Lax-Friedrichs', order=7, integrator='Runge-Kutta 4'},			-- 8.5173838468289e-15  			0.00055485267752876
	{solver='weno', wenoMethod='1996 Jiang Shu', wenoFlux='Lax-Friedrichs', order=9, integrator='Runge-Kutta 4'},			-- 7.7100869091962e-15  			0.00045701681628967
	{solver='weno', wenoMethod='2008 Borges', wenoFlux='Lax-Friedrichs', order=9, integrator='Runge-Kutta 4'},				-- 7.6095813678068e-15  			0.00050131231958936
	{solver='weno', wenoMethod='2010 Shen Zha', wenoFlux='Lax-Friedrichs', order=9, integrator='Runge-Kutta 4'},			-- 1.3609664610559e-14  			0.00051026202191686
	{solver='weno', wenoMethod='1996 Jiang Shu', wenoFlux='Lax-Friedrichs', order=11, integrator='Runge-Kutta 4'},			-- 4.4889222547417e-15  			0.00042165683328641
	{solver='weno', wenoMethod='2008 Borges', wenoFlux='Lax-Friedrichs', order=11, integrator='Runge-Kutta 4'},				-- 4.4436026039318e-15  			0.00044703551323215
	{solver='weno', wenoMethod='2010 Shen Zha', wenoFlux='Lax-Friedrichs', order=11, integrator='Runge-Kutta 4'},			-- 8.7428978987059e-15  			0.00046974951259134
	{solver='weno', wenoMethod='1996 Jiang Shu', wenoFlux='Lax-Friedrichs', order=13, integrator='Runge-Kutta 4'},			-- 7.3347361170817e-15  			0.00044939685692017
	{solver='weno', wenoMethod='2008 Borges', wenoFlux='Lax-Friedrichs', order=13, integrator='Runge-Kutta 4'},				-- 7.2548304169695e-15  			0.00046166530014287
	{solver='weno', wenoMethod='2010 Shen Zha', wenoFlux='Lax-Friedrichs', order=13, integrator='Runge-Kutta 4'},			-- 2.3371712551401e-14  			0.00047554368514062
	{solver='weno', wenoMethod='1996 Jiang Shu', wenoFlux='Lax-Friedrichs', integrator='Runge-Kutta 4, TVD'},				-- 2.4345564626771e-13  			0.00077420691307587
	{solver='weno', wenoMethod='2008 Borges', wenoFlux='Lax-Friedrichs', integrator='Runge-Kutta 4, TVD'},					-- 4.4378238063525e-14  			0.00064159945825858
	{solver='weno', wenoMethod='2010 Shen Zha', wenoFlux='Lax-Friedrichs', integrator='Runge-Kutta 4, TVD'},				-- 4.8910636825383e-14  			0.00061262672492586
	{solver='weno', wenoMethod='1996 Jiang Shu', wenoFlux='Lax-Friedrichs', integrator='Runge-Kutta 4, non-TVD'},			-- 2.4950126600171e-13  			0.00077415102265575
	{solver='weno', wenoMethod='2008 Borges', wenoFlux='Lax-Friedrichs', integrator='Runge-Kutta 4, non-TVD'},				-- 5.0651540253743e-14  			0.00064106114903538
	{solver='weno', wenoMethod='2010 Shen Zha', wenoFlux='Lax-Friedrichs', integrator='Runge-Kutta 4, non-TVD'},			-- 5.5381697491863e-14  			0.00061303349011055
	{solver='weno', wenoMethod='1996 Jiang Shu', wenoFlux='Lax-Friedrichs', integrator='Runge-Kutta 4, 3/8ths rule'},		-- 0.012520605169068    			0.0013252024019549
	{solver='weno', wenoMethod='2008 Borges', wenoFlux='Lax-Friedrichs', integrator='Runge-Kutta 4, 3/8ths rule'},			-- 0.14741644988061     			0.0040183089308677
	{solver='weno', wenoMethod='2010 Shen Zha', wenoFlux='Lax-Friedrichs', integrator='Runge-Kutta 4, 3/8ths rule'},		-- 0.99664040188526     			0.001305744394872
	{solver='weno', wenoMethod='1996 Jiang Shu', wenoFlux='Lax-Friedrichs', integrator='forward Euler'},					-- 0.023045797210864    			0.0032148202429377
	{solver='weno', wenoMethod='2008 Borges', wenoFlux='Lax-Friedrichs', integrator='forward Euler'},						-- 0.18480201977208     			0.0056315981827047
	{solver='weno', wenoMethod='2010 Shen Zha', wenoFlux='Lax-Friedrichs', integrator='forward Euler'},						-- 0.10436555718775     			0.0052821632535482
	{solver='hll', integrator='forward Euler'},																				-- 0.00037540541165354  			0.0029204942118918
	{solver='euler-burgers', integrator='forward Euler'},																	-- 0.00029551600678445  			0.0020825509821013
	{solver='euler-hllc', hllcMethod=0, integrator='forward Euler', usePLM='plm-cons-alone'},		 	  	  	  	  	  	-- 0.00039212252597743  			0.0035041915046615
	{solver='euler-hllc', hllcMethod=0, integrator='forward Euler', usePLM='plm-prim-alone'},		 	 	 	 	 	 	-- 0.0003914811986822   			0.0033715333876716
	{solver='euler-hllc', hllcMethod=0, integrator='forward Euler'},														-- 0.00029551600678424  			0.0026034369564245
	{solver='euler-hllc', hllcMethod=0, integrator='forward Euler', usePLM='plm-cons'},		 	 	 	 	 	 	 		-- 0.00029551600678424  			0.0026034369564245
	{solver='euler-hllc', hllcMethod=0, integrator='Runge-Kutta 3, TVD', usePLM='plm-prim-alone'},							-- 0.00039567583187226  			0.0038904131956502
	{solver='euler-hllc', hllcMethod=0, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons-alone'},							-- 0.00039567583187237  			0.0039987798334505
	{solver='euler-hllc', hllcMethod=0, integrator='Runge-Kutta 3, TVD'},													-- 0.00039208124005629  			0.0031547573035925
	{solver='euler-hllc', hllcMethod=0, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons'},								-- 0.00039208124005629  			0.0031547573035925
	{solver='euler-hllc', hllcMethod=1, integrator='forward Euler', usePLM='plm-cons-alone'},								-- 0.0003925245324907   			0.0035041915046615
	{solver='euler-hllc', hllcMethod=1, integrator='forward Euler', usePLM='plm-prim-alone'},								-- 0.00039148119870067  			0.0033715333876716
	{solver='euler-hllc', hllcMethod=1, integrator='forward Euler'},														-- 0.0002955160067843   			0.0026034369564245
	{solver='euler-hllc', hllcMethod=1, integrator='forward Euler', usePLM='plm-cons'},										-- 0.0002955160067843   			0.0026034369564245
	{solver='euler-hllc', hllcMethod=1, integrator='Runge-Kutta 3, TVD', usePLM='plm-prim-alone'},							-- 0.0003956758318725   			0.0038904131956501
	{solver='euler-hllc', hllcMethod=1, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons-alone'},							-- 0.00039567583187227  			0.0039987798334505
	{solver='euler-hllc', hllcMethod=1, integrator='Runge-Kutta 3, TVD'},													-- 0.00039208124005635  			0.0031547573035925
	{solver='euler-hllc', hllcMethod=1, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons'},								-- 0.00039208124005635  			0.0031547573035925
	{solver='euler-hllc', hllcMethod=2, integrator='forward Euler', usePLM='plm-cons-alone'},	 							-- 0.0003927495532771   			0.0035041915046615
	{solver='euler-hllc', hllcMethod=2, integrator='forward Euler', usePLM='plm-prim-alone'},			 	  	  	  	  	-- 0.00039148119875762  			0.0033715333876716
	{solver='euler-hllc', hllcMethod=2, integrator='forward Euler'},														-- 0.00029551600678437  			0.0026034369564245
	{solver='euler-hllc', hllcMethod=2, integrator='forward Euler', usePLM='plm-cons'},			  	  	  	  	  	  		-- 0.00029551600678437  			0.0026034369564245
	{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD', usePLM='plm-prim-alone'},	  	 	 	 	 	 	-- 0.00039567583187248  			0.0038904131956501
	{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons-alone'},							-- 0.00039567583187241  			0.0039987798334505
	{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD'},													-- 0.00039208124005633  			0.0031547573035925
	{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons'},								-- 0.00039208124005633  			0.0031547573035925
	{solver='roe', integrator='forward Euler', fluxLimiter='smart'},														-- 1.0911857176237e-06  			0.0006129313415327
	{solver='roe', integrator='forward Euler', fluxLimiter='ospre'},														-- 0.00019705480525943  			0.0020319762206545
	{solver='roe', integrator='forward Euler', fluxLimiter='Fromm'},														-- 1.5346799701183e-07  			0.00062670027933916
	{solver='roe', integrator='forward Euler', fluxLimiter='CHARM'},														-- 1.4822461340873e-06  			0.00060082731534463
	{solver='roe', integrator='forward Euler', fluxLimiter='van Albada 1'},													-- 1.2462029419466e-06  			0.00065025946311501
	{solver='roe', integrator='forward Euler', fluxLimiter='Barth-Jespersen'},												-- 4.808764902909e-07   			0.00045706811880196
	{solver='roe', integrator='forward Euler', fluxLimiter='Beam-Warming'},													-- 1.0595112307524e-06  			0.0013601751183648
	{solver='roe', integrator='forward Euler', fluxLimiter='van Leer'},														-- 7.8024435458499e-07  			0.00054262570046035
	{solver='roe', integrator='forward Euler', fluxLimiter='van Albada 2'},													-- 2.3724738671197e-06  			0.00081184935099592
	{solver='roe', integrator='forward Euler', fluxLimiter='donor cell'},													-- 0.00029551600678436  			0.0025480145819915
	{solver='roe', integrator='forward Euler', fluxLimiter='Oshker'},														-- 2.4949277440571e-06  			0.00067671846961815
	{solver='roe', integrator='forward Euler', fluxLimiter='superbee'},														-- 1.8741138235795e-06  			0.00029845606056438
	{solver='roe', integrator='forward Euler', fluxLimiter='Sweby'},														-- 2.0705889490385e-06  			0.00041546702447137
	{solver='roe', integrator='forward Euler', fluxLimiter='HQUICK'},														-- 1.3626347787182e-06  			0.00059526337064112
	{solver='roe', integrator='forward Euler', fluxLimiter='Koren'},														-- 1.0448963734558e-06  			0.00050202680619077
	{solver='roe', integrator='forward Euler', fluxLimiter='HCUS'},															-- 1.1792669058976e-06  			0.00055928895482656
	{solver='roe', integrator='forward Euler', fluxLimiter='monotized central'},											-- 4.5584837524636e-07  			0.00045225285804202
	{solver='roe', integrator='forward Euler', fluxLimiter='UMIST'},														-- 1.1751903573017e-06  			0.00061938621334895
	{solver='roe', integrator='forward Euler', fluxLimiter='minmod'},														-- 2.549668825168e-06   			0.00084753894377874
	{solver='roe', integrator='forward Euler', fluxLimiter='Lax-Wendroff'},													-- 7.5336775409067e-07  			0.0044082598083321
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='CHARM'},									-- 0.0096538616945486   			0.018033561576867
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Lax-Wendroff'},								-- 0.99794993315935     			0.22900538098487
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Beam-Warming'},								-- 0.99792029781654     			0.22642827784528
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='van Leer'},									-- 0.012965241011057    			0.24360327920558
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Fromm'},									-- 0.99786383295137     			0.22554745191931
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Barth-Jespersen'},							-- 0.99773589832872     			0.2154814116913
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='smart'},									-- 0.99593259324672     			0.095525440317095
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='superbee'},									-- 0.03588595793091     			0.0054859189942998
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='monotized central'},						-- 0.044254620253092    			0.0048356399731734
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Koren'},									-- 0.022743668518402    			0.0045319629410316
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='HQUICK'},									-- 0.01532996379064     			0.0043568816756108
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='HCUS'},										-- 0.010371614847005    			0.0034376341809745
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='donor cell'},								-- 0.00029551600678436  			0.0025480145819915
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='UMIST'},									-- 0.017416713705133    			0.0030517363635154
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Sweby'},									-- 0.0082899026003203   			0.0025329368985951
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='ospre'},									-- 0.00013263354368475  			0.0016309138954085
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Oshker'},									-- 0.0063142482387812   			0.002464024257206
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='van Albada 1'},								-- 0.0017180281251711   			0.0014345138627451
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='van Albada 2'},								-- 0.99220555463398     			0.21732912540924
	{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='minmod'},									-- 0.00053201867890155  			0.0013513761228327
	{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},													-- 0.00029592694908829  			0.0014406203418454
	{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim'},														-- 0.00029593080182222  			0.0014316307885825
	{solver='roe', integrator='forward Euler', usePLM='plm-eig'},															-- 0.00029590279585709  			0.0012223900726728
	{solver='roe', integrator='forward Euler', usePLM='plm-athena'},														-- 9.7002822784791e-05  			0.00093771140713331
	{solver='roe', integrator='forward Euler', usePLM='ppm-experimental'},													-- 0.00086335439273674  			0.24052752253483
	{solver='roe', integrator='Runge-Kutta 3', fluxLimiter='Lax-Wendroff'},													-- 0.042514748138836    			0.017088495886816
	{solver='roe', integrator='Runge-Kutta 3, TVD', fluxLimiter='Lax-Wendroff'},											-- 9.6682375641956e-05  			0.006423504743502
	{solver='roe', integrator='Runge-Kutta 4, non-TVD', fluxLimiter='Lax-Wendroff'},										-- 9.6682357228525e-05  			0.0063992773649444
	{solver='roe', integrator='Runge-Kutta 4, TVD', fluxLimiter='Lax-Wendroff'},											-- 9.6682357223262e-05  			0.0063995755451124
	{solver='roe', integrator='Runge-Kutta 4', fluxLimiter='Lax-Wendroff'},													-- 9.6682357219858e-05  			0.0063992773649444
	{solver='roe', integrator='Runge-Kutta 2, non-TVD', fluxLimiter='Lax-Wendroff'},										-- 9.6682087934282e-05  			0.27721278024955
	{solver='roe', integrator='Runge-Kutta 2, TVD', fluxLimiter='Lax-Wendroff'},											-- 9.6682087934274e-05  			0.0064919094775845
	{solver='roe', integrator='Runge-Kutta 2 Heun', fluxLimiter='Lax-Wendroff'},											-- 9.6682087934198e-05  			0.006491909477584
	{solver='roe', integrator='Runge-Kutta 2 Ralston', fluxLimiter='Lax-Wendroff'},											-- 9.6682087934305e-05  			0.0064541422444027
	{solver='roe', integrator='Runge-Kutta 2', fluxLimiter='Lax-Wendroff'},													-- 9.6682087934377e-05  			0.006417331483589
	{solver='roe', integrator='Runge-Kutta 4, 3/8ths rule', fluxLimiter='Lax-Wendroff'},									-- 2.4183176943556e-05  			0.005060892820765
	{solver='roe', integrator='Runge-Kutta 3', usePLM='plm-athena'},														-- 0.042519919280547    			0.013615227412304
	{solver='roe', integrator='Runge-Kutta 2', usePLM='plm-athena'},														-- 1.3132898617302e-06  			0.00067817159415456
	{solver='roe', integrator='Runge-Kutta 2 Ralston', usePLM='plm-athena'},												-- 1.3128886603164e-06  			0.0007200643702763
	{solver='roe', integrator='Runge-Kutta 2, TVD', usePLM='plm-athena'},													-- 1.3118684200161e-06  			0.00072699746948317
	{solver='roe', integrator='Runge-Kutta 2 Heun', usePLM='plm-athena'},													-- 1.3118684201669e-06  			0.00072699746948318
	{solver='roe', integrator='Runge-Kutta 4', usePLM='plm-athena'},														-- 1.2279744814311e-06  			0.00069783409481987
	{solver='roe', integrator='Runge-Kutta 4, non-TVD', usePLM='plm-athena'},												-- 1.2279744811973e-06  			0.0006978340948199
	{solver='roe', integrator='Runge-Kutta 4, TVD', usePLM='plm-athena'},													-- 1.2279971850421e-06  			0.00069848808747002
	{solver='roe', integrator='Runge-Kutta 3, TVD', usePLM='plm-athena'},													-- 1.226178781197e-06   			0.00070049425851808
	{solver='roe', integrator='Runge-Kutta 4, 3/8ths rule', usePLM='plm-athena'},											-- 7.2560410495091e-05  			0.00073737230229337
--[[ backward Euler is having trouble converging
	{solver='weno', wenoMethod='1996 Jiang Shu', wenoFlux='Lax-Friedrichs', order=5, integrator='backward Euler'},			-- 0.038374104574442    			0.0029568572954654
	{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=10, epsilon=1e-10}},	-- 2.6032104426088e-05  			0.0044045564570584
	{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=10, epsilon=1e-20}},	-- 1.6528741792613e-05  			0.0044050280780159
	{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=10, epsilon=1e-30}},	-- 7.4842828263317e-07  			0.0044082504092463
	{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-10}},	-- 2.6032104426088e-05  			0.0044045564570584
	{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-20}},	-- 1.0913463408278e-05  			0.0045369293796878
	{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-30}},	-- 6.7989094983795e-07  			0.0045519351868826
	{solver='roe', fluxLimiter='superbee', integrator='backward Euler', integratorArgs={restart=10, epsilon=1e-10}},		-- 2.7839818059459e-05  			0.00033572008372309
	{solver='roe', fluxLimiter='superbee', integrator='backward Euler', integratorArgs={restart=10, epsilon=1e-20}},		-- 1.9887337903574e-05  			0.0003009941921701
	{solver='roe', fluxLimiter='superbee', integrator='backward Euler', integratorArgs={restart=10, epsilon=1e-30}},		-- 1.8741498287414e-06  			0.00029845646884114
	{solver='roe', fluxLimiter='superbee', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-10}},		-- 2.7839818059459e-05  			0.00033572008372309
	{solver='roe', fluxLimiter='superbee', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-20}},		-- 8.1603905215462e-06  			0.00031412836681015
	{solver='roe', fluxLimiter='superbee', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-30}},		-- 1.8744531399e-06     			0.00030408892959842
	-- failing on advect wave
	{solver='weno', wenoMethod='2008 Borges', wenoFlux='Lax-Friedrichs', order=5, integrator='backward Euler'},				-- nan								0.0066108897724998
	{solver='weno', wenoMethod='2010 Shen Zha', wenoFlux='Lax-Friedrichs', order=5, integrator='backward Euler'},			-- nan								0.0051641313950373
--]]	
	-- failing on Sod
	--{solver='roe', integrator='Runge-Kutta 2, non-TVD', usePLM='plm-athena'},												-- 1.3161158262309e-06				nan
}
