#region Copyright
///<remarks>
/// <Graz Lagrangian Particle Dispersion Model>
/// Copyright (C) [2019]  [Dietmar Oettl, Markus Kuntner]
/// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
/// the Free Software Foundation version 3 of the License
/// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
/// You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
///</remarks>
#endregion

using System;
using System.IO;
using System.Threading;
using System.Threading.Tasks;

namespace GRAL_2001
{
    /*  GRAZ LAGRANGIAN PARTICLE MODELL GRAL
        COMPREHENSIVE DESCRIPTION CAN BE FOUND IN OETTL, 2016
        THE GRAL MODEL HAS BEEN  DEVELOPED BY DIETMAR OETTL SINCE AROUND 1999.
     */
    partial class Program
    {
        /*INPUT FILES :
          BASIC DOMAIN INFORMATION              GRAL.geb
          MAIN CONTROL PARAM FILEs              in.dat
          GEOMETRY DATA                         ggeom.asc
          LANDUSE DATA                          landuse.asc
          METEOROLOGICAL DATA                   meteopgt.all, or inputzr.dat, or sonic.dat
          MAX. NUMBER OF CPUs                   Max_Proc.txt
          EMISSION SOURCES                      line.dat, point.dat, cadastre.dat, portals.dat
          TUNNEL JET DESTRUCTION
          BY TRAFFIC ON OPPOSITE LANE           oppsite_lane.txt
          POLLUTANTS SUCKED IN BY
          TUNNEL PORTAL AT OPPOSITE LANE        tunnel_entrance.txt
          LOCATIONS AND HEIGHTS OF BUILDINGS    buildings.dat
          LOCATIONS AND HEIGHTS OF VEGETATION   vegetation.dat
          LOCATIONS AND HEIGHTS OF RECEPTORS    Receptor.dat
          NUMBER OF VERTICAL LAYERS FOR
          PROGNOSTIC FLOW FIELD MODEL           micro_vert_layers.txt
          RELAXATION FACTORS FOR
          PROGNOSTIC FLOW FIELD MODEL           relaxation_factors.txt
          MINIMUM ANd MAXIMUM INTEGRATION TIMES
          FOR PROGNOSTIC FLOW FIELD MODEL       Integrationtime.txt
          ROUGHNESS LENGTH FOR OBSTACLES
          FOR PROGNOSTIC FLOW FIELD MODEL       building_roughness.txt
          TRANSIENT GRAL MODE CONC.THRESHOLD	GRAL_Trans_Conc_Threshold.txt
          POLLUTANT & WET DEPOSITION SETTINGS	Pollutant.txt
          WET DEPOSITION PRECIPITATION DATA		Precipitation.txt
          CALCULATE 3D CONCENTRATION DATA IN
          TRANSIENT GRAL MODE					GRAL_Vert_Conc.txt
         */

        // 鏀惧湪 partial class Program 閲岋紙浠绘剰瀛楁鍖洪兘琛岋級
        private static readonly object _ioLock = new object();
        public static void WriteLineSafe(string path, string line, bool append = true)
        {
            lock (_ioLock)
            {
                // 鍏佽鍏跺畠璇汇€佸苟鍙戣拷鍔狅紝閬垮厤鍏变韩鍐茬獊
                using var fs = new FileStream(path,
                    append ? FileMode.Append : FileMode.Create,
                    FileAccess.Write,
                    FileShare.ReadWrite);
                using var sw = new StreamWriter(fs);
                sw.WriteLine(line);
                sw.Flush();
            }
        }

        public static float GetEmissionFactorForSourceGroup(int iwet, int realSourceGroup, float fallbackFactor) //20260521
        {
            if (!EmissionTimeseriesExist || EmFacTimeSeries == null)
            {
                return fallbackFactor;
            }

            int row = iwet - 1;
            int sgInternal = SourceGroups.IndexOf(realSourceGroup);
            if (row < 0 || sgInternal < 0 ||
                row > EmFacTimeSeries.GetUpperBound(0) ||
                sgInternal > EmFacTimeSeries.GetUpperBound(1))
            {
                return 0f;
            }

            return EmFacTimeSeries[row, sgInternal];
        }

        private static float GetMaxEmissionFactorForCurrentStep() //20260521
        {
            if (!EmissionTimeseriesExist || EmFacTimeSeries == null)
            {
                return 1f;
            }

            int row = IWET - 1;
            if (row < 0 || row > EmFacTimeSeries.GetUpperBound(0))
            {
                return 0f;
            }

            float maxFactor = 0f;
            int sgMax = Math.Min(SourceGroups.Count - 1, EmFacTimeSeries.GetUpperBound(1));
            for (int sg = 0; sg <= sgMax; sg++)
            {
                maxFactor = MathF.Max(maxFactor, EmFacTimeSeries[row, sg]);
            }

            return maxFactor;
        }

        public static DenseGasState CreateInitialDenseGasState() //20260521
        {
            float rhoAir0 = 1.2041f * (273.15f / MathF.Max(200f, AmbientTempK));
            float rhoGas0 = (GasDensityRatio0 > 0.0f && GasDensityRatio0 != 1.0f)
                ? rhoAir0 * GasDensityRatio0
                : rhoAir0 * (GasMolWeight / 29.0f) * (AmbientTempK / MathF.Max(150f, GasTempK));

            return new DenseGasState
            {
                rhoBulk = rhoGas0,
                cloudDepth = _initCloudDepth,
                halfWidthX = _initCloudRadius,
                halfWidthY = _initCloudRadius,
                up = 0f,
                vp = 0f,
                wp = 0f,
                phase = DensePhase.Descent,
                spreadAngle = 0f,
                Ug = 0f,
                spreadActive = false,
                hitGroundOnce = false,
                HbulkLocal = 0f,
                RhoBulkLocal = 0f,
                BbulkLocal = 0f,
                ColumnCount = 0,
                dbgColHbulk = 0f,
                dbgColRhoBulk = 0f,
                dbgColBbulk = 0f
            };
        }

        static void Main(string[] args)
        {
            int p = (int)Environment.OSVersion.Platform;
            if ((p == 4) || (p == 6) || (p == 128))
            {
                //Console.WriteLine ("Running on Unix");
                RunOnUnix = true;
            }

            //WRITE GRAL VERSION INFORMATION TO SCREEN
            Console.WriteLine("");
            Console.WriteLine("+------------------------------------------------------+");
            Console.WriteLine("|                                                      |");
            string Info = "+  > >         G R A L VERSION: 24.11            < <   +";
            Console.WriteLine(Info);
            if (RunOnUnix)
            {
                Console.WriteLine("|                     L I N U X                        |");
            }
#if NET6_0
            Console.WriteLine("|                   .NET6 Version                      |");
#elif NET7_0
            Console.WriteLine("|                   .NET7 Version                      |");
#elif NET8_0_OR_GREATER
            Console.WriteLine("|                   .NET8 Version                      |");
#else
            Console.WriteLine("|                 .Net Core Version                    |");
#endif
            Console.WriteLine("|                                                      |");
            Console.WriteLine("+------------------------------------------------------+");
            Console.WriteLine("");

            ShowCopyright(args);

            // write zipped files?
            ResultFileZipped = false;

            LogLevel = CheckCommandLineArguments(args);

            //Delete file Problemreport_GRAL.txt
            if (File.Exists("Problemreport_GRAL.txt") == true)
            {
                try
                {
                    File.Delete("Problemreport_GRAL.txt");
                }
                catch { }
            }
            // Write to "Logfile_GRALCore"
            try
            {
                ProgramWriters.LogfileGralCoreWrite(new String('-', 80));
                ProgramWriters.LogfileGralCoreWrite(Info);
                ProgramWriters.LogfileGralCoreWrite("Computation started at: " + DateTime.Now.ToString());
                ProgramWriters.LogfileGralCoreWrite("Computation folder:     " + Directory.GetCurrentDirectory());
                Info = "Application hash code:  " + GetAppHashCode();
                ProgramWriters.LogfileGralCoreWrite(Info);
            }
            catch { }

            ProgramReaders ReaderClass = new ProgramReaders();

            //Read number of user defined vertical layers
            ReaderClass.ReadMicroVertLayers();
            GFFFilePath = ReaderClass.ReadGFFFilePath();

            //Lowest grid level of the prognostic flow field model grid
            HOKART[0] = 0;

            //read main GRAL domain file "GRAL.geb" and check if the file "UseOrigTopography.txt" exists -> needed to read In.Dat!
            ReaderClass.ReadGRALGeb();
            Console.WriteLine($"[DBG] after in.dat XsiMin={Program.XsiMinGral} XsiMax={Program.XsiMaxGral} EtaMin={Program.EtaMinGral} EtaMax={Program.EtaMaxGral}");//20260201

            //optional: read GRAMM file ggeom.asc
            if (File.Exists("ggeom.asc") == true)
            {
                ReaderClass.ReadGRAMMGeb();
            }

            //horizontal grid sizes of the GRAL concentration grid
            GralDx = (float)((XsiMaxGral - XsiMinGral) / (float)NXL);
            GralDy = (float)((EtaMaxGral - EtaMinGral) / (float)NYL);

            //Set the maximum number of threads to be used in each parallelized region
            ReaderClass.ReadMaxNumbProc();

            //sets the number of cells near the walls of obstacles, where a boundary-layer is computed in the diagnostic flow field approach
            IGEB = Math.Max((int)(20 / DXK), 1);

            //Read Pollutant and Wet deposition data
            Odour = ReaderClass.ReadPollutantTXT();

            // === Read dense/light gas overrides (A-scheme) ===
            ReadDenseLightConfig("densegas.cfg");

            // === Build fingerprint (write to run directory) ===
            try
            {
                var asm = System.Reflection.Assembly.GetExecutingAssembly().Location;
                string cwd = Environment.CurrentDirectory;
                string fp  = System.IO.Path.Combine(cwd, "gral_build_fingerprint.txt");
                System.IO.File.WriteAllLines(fp, new []
                {
                    "=== GRAL Dense/Light Build Fingerprint ===",
                    $"Time={DateTime.Now:yyyy-MM-dd HH:mm:ss}",
                    $"Assembly={asm}",
                    $"GasBuoyancyMode={GasBuoyancyMode}",
                    $"UseDenseGasMicroSpray={UseDenseGasMicroSpray}",
                    $"AmbientTempK={AmbientTempK}",
                    $"GasTempK={GasTempK}",
                    $"GasDensityRatio0={GasDensityRatio0}",
                    $"GasMolWeight={GasMolWeight}",
                    $"DenseMixFactor={DenseMixFactor}",
                    $"DenseEntrainmentFactor={DenseEntrainmentFactor}",
                    $"UseExplicitExitVelocityDirection={UseExplicitExitVelocityDirection}" //20260521 1852 log Trial 8 explicit point-source velocity switch
                });
            }
            catch { /* no-op */ }
      
            Console.WriteLine("=== GRAL Dense/Light Build ===");
            Console.WriteLine($"Assembly: {System.Reflection.Assembly.GetExecutingAssembly().Location}");
            Console.WriteLine($"GasBuoyancyMode={GasBuoyancyMode}, UseMicroSpray={UseDenseGasMicroSpray}");
            Console.WriteLine($"GasDensityRatio0={GasDensityRatio0}, DenseVertCoeff={DenseVertCoeff}, DenseMixFactor={DenseMixFactor}, DenseEntrainmentFactor={DenseEntrainmentFactor}");
            Console.WriteLine($"UseExplicitExitVelocityDirection={UseExplicitExitVelocityDirection}"); //20260521 1852 expose explicit point-source velocity switch at startup

            //Create large arrays
            CreateLargeArrays();

            //optional: reading GRAMM orography file ggeom.asc -> there are two ways to read ggeom.asc (1) the files contains all information or (2) the file just provides the path to the original ggeom.asc
            ReaderClass.ReadGgeomAsc();

            //number of grid cells of the GRAL microscale flow field
            NII = (int)((XsiMaxGral - XsiMinGral) / DXK);
            NJJ = (int)((EtaMaxGral - EtaMinGral) / DYK);
            AHKOri = CreateArray<float[]>(NII + 2, () => new float[NJJ + 2]);
            GralTopofile = ReaderClass.ReadGRALTopography(NII, NJJ); // GRAL Topofile OK?

            //reading main control file in.dat
            ReaderClass.ReadInDat();
            //total number of particles released for each weather situation
            NTEILMAX = (int)(TAUS * TPS);
            //Volume of the GRAL concentration grid
            GridVolume = GralDx * GralDy * GralDz;

            //Reading building data
            //case 1: complex terrain
            {
                int SIMD = CheckSIMD();
                if (Topo == Consts.TerrainAvailable)
                {
                    InitGralTopography(SIMD);
                    ReaderClass.ReadBuildingsTerrain(Program.CUTK); //define buildings in GRAL
                }
                //flat terrain application
                else
                {
                    InitGralFlat();
                    ReaderClass.ReadBuildingsFlat(Program.CUTK); //define buildings in GRAL
                }
            }

            // array declarations for prognostic and diagnostic flow field
            if ((FlowFieldLevel > Consts.FlowFieldNoBuildings) || (Topo == Consts.TerrainAvailable))
            {
                // create jagged arrays manually to keep memory areas of similar indices togehter -> reduce false sharing & 
                // save memory because of the unused index 0  
                DIV = new float[NII + 1][][];
                DPM = new float[NII + 1][][];
                for (int i = 1; i < NII + 1; ++i)
                {
                    DIV[i] = new float[NJJ + 1][];
                    DPM[i] = new float[NJJ + 1][];
                    for (int j = 1; j < NJJ + 1; ++j)
                    {
                        DIV[i][j] = new float[NKK + 1];
                        DPM[i][j] = new float[NKK + 2];
                    }
                    if (i % 100 == 0)
                    {
                        Console.Write(".");
                    }
                }
            }

            //In case of transient simulations: load presets and define arrays
            if (ISTATIONAER == Consts.TransientMode)
            {
                TransientPresets.LoadAndDefine();
                ReaderClass.Read3DSnapshotSettings(); // optional: GRAL_3D_Snapshots.txt    20260401add
            }

            Console.WriteLine(".");
            //reading receptors from file Receptor.dat
            ReaderClass.ReadReceptors();

            //the horizontal standard deviations of wind component fluctuations are dependent on the averaging time (dispersion time)
            if ((IStatistics == Consts.MeteoPgtAll))
            {
                StdDeviationV = (float)Math.Pow(TAUS / 3600, 0.2);
            }

            //for applications in flat terrain, the roughness length is homogenous as defined in the file in.dat
            if (Topo != Consts.TerrainAvailable)
            {
                Z0Gramm[1][1] = Z0;
            }

            //checking source files
            if (File.Exists("line.dat") == true)
            {
                LS_Count = 1;
            }

            if (File.Exists("portals.dat") == true)
            {
                TS_Count = 1;
            }

            if (File.Exists("point.dat") == true)
            {
                PS_Count = 1;
            }

            if (File.Exists("cadastre.dat") == true)
            {
                AS_Count = 1;
            }
            
            if (File.Exists("volume.dat") == true)//20260310
            {
                VS_Count = 1;
            }

            Console.WriteLine();
            Info = "Total number of horizontal slices for concentration grid: " + NS.ToString();
            Console.WriteLine(Info);

            for (int i = 0; i < NS; i++)
            {
                try
                {
                    Info = "  Slice height above ground [m]: " + HorSlices[i].ToString();
                    Console.WriteLine(Info);
                }
                catch
                { }
            }

            //emission modulation for transient mode
            ReaderClass.ReadEmissionTimeseries();

            //Create Vegetation array only in areas of interest and at point 0/0 -> reduce memory footprint
            //VEG = CreateArray<float[][]>(NII + 2, () => CreateArray<float[]>(NJJ + 2, () => new float[NKK + 1]));
            VEG = new float[NII + 2][][];
            for (int i = 0; i < NII + 2; i++)
            {
                VEG[i] = new float[NJJ + 2][];
                if (i == 0)
                {
                    VEG[i][0] = new float[NKK + 1];
                }
            }
            Program.VEG[0][0][0] = -0.001f; //sign for reading vegetation one times, when claculating flow fields
            Program.VEG[0][0][1] = -0.001f; //sign for reading vegetation, if *.gff files exist
            COV = CreateArray<float[]>(NII + 2, () => new float[NJJ + 2]);

            //reading optional files used to define areas where either the tunnel jet stream is destroyed due to traffic
            //on the opposite lanes of a highway or where pollutants are sucked into a tunnel portal
            ReaderClass.ReadTunnelfilesOptional();

            //optional: reading land-use data from file landuse.asc
            ReaderClass.ReadLanduseFile();

            //Use the adaptive roughness lenght mode?
            if (AdaptiveRoughnessMax > 0 && BuildingsExist)
            {
                //define FlowField dependend Z0, UStar and OL, generate Z0GRAL[][] array and write RoghnessGRAL file
                InitAdaptiveRoughnessLenght(ReaderClass);
            }
            else
            {
                AdaptiveRoughnessMax = 0;
            }

            if (Program.UseFixedRndSeedVal)
            {
                Info = "GRAL Deterministic Mode";
                Console.WriteLine("[PG] reached MODE branch Program.cs");
                Console.WriteLine(Info);
                ProgramWriters.LogfileGralCoreWrite(Info);
            }

            //setting the lateral borders of the domain -> Attention: changes the GRAL borders from absolute to relative values!
            InitGralBorders();
            Console.WriteLine($"[DBG] after InitGralBorders XsiMin={Program.XsiMinGral} XsiMax={Program.XsiMaxGral} EtaMin={Program.EtaMinGral} EtaMax={Program.EtaMaxGral}");//20260202

            //reading point source data
            if (PS_Count > 0)
            {
                PS_Count = 0;
                Console.WriteLine();
                Console.WriteLine("Reading file point.dat");
                ReadPointSources.Read();
            }
            //reading line source data
            if (LS_Count > 0)
            {
                LS_Count = 0;
                Console.WriteLine();
                Console.WriteLine("Reading file line.dat");
                ReadLineSources.Read();
            }
            //reading tunnel portal data
            if (TS_Count > 0)
            {
                TS_Count = 0;
                Console.WriteLine();
                Console.WriteLine("Reading file portals.dat");
                ReadTunnelPortals.Read();
            }
            //reading area source data
            if (AS_Count > 0)
            {
                AS_Count = 0;
                Console.WriteLine();
                Console.WriteLine("Reading file cadastre.dat");
                ReadAreaSources.Read();
            }
            //reading volume source data//20260310
            if (VS_Count > 0)
            {
                VS_Count = 0;
                Console.WriteLine();
                Console.WriteLine("Reading file volume.dat");
                ReadVolumeSources.Read();
                if (Program.VS_Count > 0)
                {
                    double qCont = 0.0;
                    double qInst = 0.0;
                    for (int iv = 1; iv <= Program.VS_Count; iv++)
                    {
                        if (Program.VS_Instant[iv]) qInst += Program.VS_TotalRelease[iv];
                        else qCont += Program.VS_ER[iv];
                    }
                    Console.WriteLine($"[VS] count={Program.VS_Count}, continuous={qCont:0.###} kg/h, instant_total={qInst:0.###} kg");
                }
            }

            Console.WriteLine($"[DBG] ContinuousTraj={Program.ContinuousTraj}");//20260202dbg

            //ReaderClass.RemovePrognosticSubDomainsFarDistanceToSources();//20260201娉ㄩ噴
            if (!Program.ContinuousTraj)
            {
                Console.WriteLine("[DBG] shrink subdomain: ON");//20260202
                ReaderClass.RemovePrognosticSubDomainsFarDistanceToSources();
            }
            else//20260202
            {
                Console.WriteLine("[DBG] shrink subdomain: OFF");
            }
            Console.WriteLine($"[DBG] after RemovePrognosticSubDomainsFarDistanceToSources XsiMax={Program.XsiMaxGral} EtaMax={Program.EtaMaxGral}");//20260201


            // Read transient meteorology before particle capacity is fixed. //20260521
            InputMettimeSeries InputMetTimeSeries = new InputMettimeSeries(); //20260521
            if ((ISTATIONAER == Consts.TransientMode) && (IStatistics == Consts.MeteoPgtAll)) //20260521
            {
                InputMetTimeSeries.WindData = MeteoTimeSer;
                InputMetTimeSeries.ReadMetTimeSeries();
                MeteoTimeSer = InputMetTimeSeries.WindData;
                if (MeteoTimeSer.Count == 0) // no data available
                {
                    string err = "Error when reading file mettimeseries.dat -> Execution stopped: press ESC to stop";
                    Console.WriteLine(err);
                    ProgramWriters.LogfileProblemreportWrite(err);
                }
                InputMetTimeSeries.WindData = MeteopgtLst;
                InputMetTimeSeries.ReadMeteopgtAll();
                MeteopgtLst = InputMetTimeSeries.WindData;
                if (MeteopgtLst.Count == 0) // no data available
                {
                    string err = "Error when reading file meteopgt.all -> Execution stopped: press ESC to stop";
                    Console.WriteLine(err);
                    ProgramWriters.LogfileProblemreportWrite(err);
                }
                // Read precipitation data
                ReaderClass.ReadPrecipitationTXT();
            }
            else
            {
                WetDeposition = false;
            }

            TransientReleaseEnabled = (ISTATIONAER == Consts.TransientMode) &&
                                      ContinuousTraj &&
                                      EmissionTimeseriesExist &&
                                      MeteoTimeSer.Count > 0; //20260521

            //  A1
            Console.WriteLine("[PG] before ParticleManagement.Calculate()");
            //distribution of all particles over all sources according to their source strengths (the higher the emission rate the larger the number of particles)
            NTEILMAX = ParticleManagement.Calculate();
            Console.WriteLine($"[PG] after ParticleManagement.Calculate(): NTEILMAX={NTEILMAX}");
            if (NTEILMAX < 0)
            {
                Console.WriteLine("[PG] assignment failed -> abort");
                return; 
            }


            if (GasBuoyancyMode != 0)//20251222 1736
            {
                // === INIT Dense/Light gas bulk states ===
            DenseStates = new DenseGasState[NTEILMAX + 2]; //20260521
            Console.WriteLine($"[PG] init DenseStates NTEILMAX={NTEILMAX}");
            for (int ip = 1; ip <= NTEILMAX; ip++)
            {
                DenseStates[ip] = CreateInitialDenseGasState(); //20260521
                }
            DenseStatesTransientBySG = new DenseGasState[Math.Max(1, SourceGroups.Count + 1)];
            for (int isg = 0; isg < DenseStatesTransientBySG.Length; isg++)
            {
                DenseStatesTransientBySG[isg] = CreateInitialDenseGasState(); //20260521
            }
            InitializeDenseColumnFields();//20260307
            Console.WriteLine("[PG] init DenseStates done");
            }

            //Coriolis parameter
            CorolisParam = (float)(2 * 7.29 * 0.00001 * Math.Sin(Math.Abs(LatitudeDomain) * 3.1415 / 180));

            //weather situation to start with
            IWET = IWETstart - 1;

            // freeze time base once; do not follow later IWETstart rewrites//20260402
            SnapshotIWETBase = IWETstart; 

            // read vegetation deposition factors
            VegetationDepoVelFactors = ReaderClass.ReadVegetationDeposition();

            if (ISTATIONAER == Consts.TransientMode)
            {
                Info = "Transient GRAL mode. Number of weather situations: " + MeteoTimeSer.Count.ToString();
                Console.WriteLine("[PG] reached MODE branch Program.cs");
                Console.WriteLine(Info);
                ProgramWriters.LogfileGralCoreWrite(Info);
            }

            /********************************************
             * MAIN LOOP OVER EACH WEATHER SITUATION    *
             ********************************************/
            Thread ThreadWriteGffFiles = null;
            Thread ThreadWriteConz4dFile = null;
            Thread ThreadWrite2DConcentrationFiles = null;
            int recentWeatherSituation = 0;

            bool FirstLoop = true;
            ActiveParticleMax = 0; //20260521
            while (IEND == Consts.CalculationRunning)
            {
                // if GFF writing thread has been started -> wait until GFF--WriteThread has been finished
                if (ThreadWriteGffFiles != null)
                {
                    ThreadWriteGffFiles.Join(5000); // wait up to 5s or until thread has been finished
                    if (ThreadWriteGffFiles.IsAlive)
                    {
                        Console.Write("Writing *.gff file..");
                        while (ThreadWriteGffFiles.IsAlive)
                        {
                            ThreadWriteGffFiles.Join(30000); // wait up to 30s or until thread has been finished
                            Console.Write(".");
                        }
                        Console.WriteLine();
                    }
                    ThreadWriteGffFiles = null; // Release ressources				
                }

                //Next weather situation
                IWET++;
                IDISP = IWET;
                WindVelGramm = -99; WindDirGramm = -99; StabClassGramm = -99; // Reset GRAMM values

                //show actual computed weather situation
                Console.WriteLine("_".PadLeft(79, '_'));
                Console.Write("Weather number: " + IWET.ToString());
                if (ISTATIONAER == Consts.TransientMode)
                {
                    int _iwet = IWET - 1;
                    if (_iwet < MeteoTimeSer.Count)
                    {
                        //Console.WriteLine(" - " + MeteoTimeSer[_iwet].Day + "." + MeteoTimeSer[_iwet].Month + "-" + MeteoTimeSer[_iwet].Hour + ":00");
                        var wd = MeteoTimeSer[_iwet];//20260101 HH:MM:SS
                        string timeStr = $"{wd.HourInt:D2}:{wd.Minute:D2}:{wd.Second:D2}";
                        Console.WriteLine(" - " + wd.Day + "." + wd.Month + "-" + timeStr);
                    }
                }
                else
                {
                    Console.WriteLine();
                }

                //time stamp to evaluate computation times
                int StartTime = Environment.TickCount;

                //Set the maximum number of threads to be used in each parallelized region
                ReaderClass.ReadMaxNumbProc();

                //read meteorological input data 
                if (IStatistics == Consts.MeteoPgtAll)
                {
                    if (ISTATIONAER == Consts.TransientMode)
                    {
                        //search for the corresponding weather situation in meteopgt.all
                        IDISP = InputMetTimeSeries.SearchWeatherSituation() + 1;
                        IWETstart = IDISP;
                        if (IWET > MeteoTimeSer.Count)
                        {
                            IEND = Consts.CalculationFinished;
                        }
                    }
                    else
                    {
                        //in stationary mode, meteopgt.all is read here, because we need IEND
                        IEND = Input_MeteopgtAll.Read();
                        IWETstart--; // reset the counter, because metepgt.all is finally read in ReadMeteoData, after the GRAMM wind field is available
                    }
                }

                if (ISTATIONAER == Consts.TransientMode && IDISP == 0)
                {
                    break; // reached last line in mettimeseries -> exit loop				
                }

                if (IEND == Consts.CalculationFinished)
                {
                    break; // reached last line in meteopgt.all ->  exit loop
                }

                if (Program.IWET == 1)//20260201
                {
                    Console.WriteLine($"[DBG] before dispersion XsiMin={Program.XsiMinGral} XsiMax={Program.XsiMaxGral} EtaMin={Program.EtaMinGral} EtaMax={Program.EtaMaxGral}");
                }

                String WindfieldPath = ReadWindfeldTXT();

                if (Topo == Consts.TerrainAvailable)
                {
                    Console.Write("Reading GRAMM wind field: ");
                }

                if (ISTATIONAER == Consts.TransientMode && IDISP < 0) // found no corresponding weather situation in meteopgt.all
                {
                    Info = "Cannot find a corresponding entry in meteopgt.all to the following line in mettimeseries.dat - situation skipped: " + IWET.ToString();
                    Console.WriteLine();
                    Console.WriteLine(Info);
                    ProgramWriters.LogfileGralCoreWrite(Info);
                }
                else if (Topo == Consts.TerrainFlat || (Topo == Consts.TerrainAvailable && ReadWndFile.Read(WindfieldPath))) // stationary mode or an entry in meteopgt.all exist && if topo -> wind field does exist
                {
                    //GUI output
                    try
                    {
                        using (StreamWriter wr = new StreamWriter("DispNr.txt"))
                        {
                            wr.WriteLine(IWET.ToString());
                        }
                    }
                    catch { }

                    //Topography mode -> read GRAMM stability classes
                    if (Topo == Consts.TerrainAvailable)
                    {
                        string SclFileName = String.Empty;
                        if (WindfieldPath == String.Empty)
                        {
                            //case 1: stability fields are located in the same sub-directory as the GRAL executable
                            SclFileName = Convert.ToString(Program.IDISP).PadLeft(5, '0') + ".scl";
                        }
                        else
                        {
                            //case 2: wind fields are imported from a different project
                            SclFileName = Path.Combine(WindfieldPath, Convert.ToString(Program.IDISP).PadLeft(5, '0') + ".scl");
                        }

                        if (File.Exists(SclFileName))
                        {
                            Console.WriteLine("Reading GRAMM stability classes: " + SclFileName);

                            ReadSclUstOblClasses Reader = new ReadSclUstOblClasses
                            {
                                FileName = SclFileName,
                                Stabclasses = AKL_GRAMM
                            };
                            Reader.ReadSclFile();
                            AKL_GRAMM = Reader.Stabclasses;
                            Reader.Close();
                        }
                    }

                    //Read meteorological input data depending on the input type IStatisitcs
                    ReadMeteoData();

                    //Check if all weather situations have been computed
                    if (IEND == Consts.CalculationFinished)
                    {
                        break;
                    }

                    //potential temperature gradient
                    CalculateTemperatureGradient();

                    //optional: read pre-computed GRAL flow fields
                    bool GffFiles = ReadGralFlowFields.Read();

                    if (GffFiles == false) // Reading of GRAL Flowfields not successful
                    {
                        //microscale flow field: complex terrain
                        if (Topo == Consts.TerrainAvailable)
                        {
                            MicroscaleTerrain.Calculate(FirstLoop);
                        }
                        //microscale flow field: flat terrain
                        else if ((Topo == Consts.TerrainFlat) && ((BuildingsExist == true) || (File.Exists("vegetation.dat") == true)))
                        {
                            MicroscaleFlat.Calculate();
                        }
                        if (Topo == Consts.TerrainFlat)
                        {
                            Program.AHMIN = 0;
                        }

                        //optional: write GRAL flow fields
                        ThreadWriteGffFiles = new Thread(WriteGRALFlowFields.Write);
                        ThreadWriteGffFiles.Start(); // start writing thread

                        if (FirstLoop)
                        {
                            WriteGRALFlowFields.WriteGRALGeometries();
                        }
                    }
                    else
                    {
                        //Read Vegetation areas
                        if ((Program.FlowFieldLevel == Consts.FlowFieldProg) && (File.Exists("vegetation.dat") == true) && Program.VEG[0][0][1] < 0)
                        {
                            //read vegetation only once
                            ProgramReaders Readclass = new ProgramReaders();
                            Readclass.ReadVegetation();
                            Program.VEG[0][0][1] = 0;
                            Program.VEG[0][0][0] = -0.001f; //reset marker, that vegetation is available
                        }
                    }
                    ProgramWriters WriteClass = new ProgramWriters();

                    //time needed for wind-field computations
                    double CalcTimeWindfield = (Environment.TickCount - StartTime) * 0.001;

                    //reset receptor concentrations
                    if (ReceptorsAvailable) // if receptors are acitvated
                    {
                        ReadReceptors.ReceptorResetConcentration();
                    }

                    // if writing thread for *.grz files has been started -> wait until Thread has been finished
                    if (ThreadWrite2DConcentrationFiles != null)
                    {
                        ThreadWrite2DConcentrationFiles.Join(5000); // wait up to 5s or until thread has been finished
                        if (ThreadWrite2DConcentrationFiles.IsAlive)
                        {
                            Console.Write("Finish writing *.grz file..");
                            while (ThreadWrite2DConcentrationFiles.IsAlive)
                            {
                                ThreadWrite2DConcentrationFiles.Join(30000); // wait up to 30s or until thread has been finished
                                Console.Write(".");
                            }
                            Console.WriteLine();
                        }
                        ThreadWrite2DConcentrationFiles = null; // Release ressources				
                    }

                    if (FirstLoop)
                    {
                        //read receptors
                        ReceptorNumber = 0;
                        if (ReceptorsAvailable) // if receptors are acitvated
                        {
                            ReadReceptors.ReadReceptor(); // read coordinates of receptors - flow field data needed
                        }
                        //in case of complex terrain and/or the presence of buildings some data is written for usage in the GUI (visualization of vertical slices)
                        WriteClass.WriteGRALGeometries();
                        //optional: write building heights as utilized in GRAL
                        WriteClass.WriteBuildingHeights("building_heights.txt", Program.BUI_HEIGHT, "0.0", 1, Program.IKOOAGRAL, Program.JKOOAGRAL);
                        //optional: write sub Domains as utilized in GRAL
                        WriteClass.WriteSubDomain("PrognosticSubDomainAreas.txt", Program.ADVDOM, "0", 1, Program.IKOOAGRAL, Program.JKOOAGRAL);
                    }

                    RnGSeed = new DeterministicRandomGenerator(IWET, WindVelGral, WindDirGral);

                    //calculating momentum and bouyancy forces for point sources
                    PointSourceHeight.CalculatePointSourceHeight();
/*
                    //defining the initial positions of particles
                    StartCoordinates.Calculate();
*///20260113
                    int releasedThisStep = 0; //20260521
                    float releaseEmissionFactor = 1f; //20260521
                    if (Program.TransientReleaseEnabled) //20260521
                    {
                        releaseEmissionFactor = GetMaxEmissionFactorForCurrentStep(); //20260521
                        if (releaseEmissionFactor > 0f)
                        {
                            int releaseStart = ActiveParticleMax + 1;
                            int releaseEnd = releaseStart + ParticlesPerReleaseStep;
                            if (releaseEnd - 1 > TotalParticleCapacity)
                            {
                                string err = $"Transient continuous release capacity exceeded at IWET={IWET}: requested={releaseEnd - 1}, capacity={TotalParticleCapacity}";
                                Console.WriteLine(err);
                                ProgramWriters.LogfileProblemreportWrite(err);
                                IEND = Consts.CalculationFinished;
                                break;
                            }

                            StartCoordinates.Calculate(releaseStart, releaseEnd, ParticlesPerReleaseStep, releaseEmissionFactor);
                            ActiveParticleMax = releaseEnd - 1;
                            releasedThisStep = releaseEnd - releaseStart;
                        }
                    }
                    else if (Program.ContinuousTraj)//20260113 //20260521
                    {
                        if (FirstLoop)
                        {
                            StartCoordinates.Calculate();
                            ActiveParticleMax = NTEILMAX; //20260521
                        }
                    }
                    else
                    {
                        StartCoordinates.Calculate();
                        ActiveParticleMax = NTEILMAX; //20260521
                    }

                    //boundary-layer height
                    CalculateBoudaryLayerHeight();

                    if (GasBuoyancyMode != 0)
                    {
                        ResetDenseColumnFields();
                    }

                    //show meteorological parameters on screen
                    OutputOfMeteoData();
                    Console.WriteLine($"[DBG] FF grid: NII={Program.NII} NJJ={Program.NJJ} IKOOAGRAL={Program.IKOOAGRAL} JKOOAGRAL={Program.JKOOAGRAL} DXK={Program.DXK} DYK={Program.DYK}");
                    
                    //Transient mode: non-steady-state particles
                    if (ISTATIONAER == Consts.TransientMode)
                    {
                        Console.WriteLine();
                        Console.Write("Dispersion computation.....");
                        Console.WriteLine("[PG#1] before DISPERSION Program.cs"); 
                        // calculate wet deposition parameter
                        if (IWET < WetDepoPrecipLst.Count)
                        {
                            if (WetDepoPrecipLst[IWET] > 0)
                            {
                                WetDepoRW = Wet_Depo_CW * Math.Pow(WetDepoPrecipLst[IWET], WedDepoAlphaW);
                                WetDepoRW = Math.Max(0, Math.Min(1, WetDepoRW));
                            }
                            else
                            {
                                WetDepoRW = 0;
                            }
                        }

                        //set lower concentration threshold for memory effect
                        TransConcThreshold = ReaderClass.ReadTransientThreshold();
                        int cellNr = NII * NJJ;
                        DispTimeSum = TAUS;
                        // start the calculation for transient particle from all cells
                        //ParallelTransientParticleDriver(0, cellNr);20260113
                        if (!Program.ContinuousTraj)//20260113
                        {
                            ParallelTransientParticleDriver(0, cellNr);
                        }
                    } // non-steady-state particles

                    Console.WriteLine();
                    Console.Write("Dispersion computation.....");
                    Console.WriteLine("[PG#2] before DISPERSION Program.cs");  //new released the particles from all sources
                    DispTimeSum = TAUS;
                    
                    if (Program.ContinuousTraj && IWET % 10 == 0)//20260114
                    {
                        int alive = 0;
                        int activeMaxForLog = Program.TransientReleaseEnabled ? Program.ActiveParticleMax : NTEILMAX; //20260521
                        for (int i = 1; i <= activeMaxForLog; i++) //20260521
                        {
                            if (Program.ParticleSource[i] != 0 && Program.ParticleMass[i] > 0)
                                alive++;
                        }
                        Console.WriteLine($"[DBG] IWET={IWET} alive={alive}");
                        Console.WriteLine($"[DBG] rm out={Program.Remove_OutDomain} ground={Program.Remove_Ground} top={Program.Remove_Top} mass0={Program.Remove_MassZero}");//20260114
                        Console.WriteLine($"[DBG] rm other={Program.Remove_Other} keep={Program.Keep_EndOfStep} ffcell={Program.Remove_FFCell} tunnel={Program.Remove_Tunnel} maxloops={Program.Remove_MaxLoops}");//20260114
                    }
                    
                    int particleUpperBound = Program.TransientReleaseEnabled ? ActiveParticleMax + 1 : NTEILMAX + 1; //20260521
                    if (Program.TransientReleaseEnabled) //20260521
                    {
                        if (particleUpperBound > 1)
                        {
                            ParallelParticleDriver(1, particleUpperBound); //20260521
                        }
                    }
                    else if (Program.ContinuousTraj)//20260113 //20260521
                    {
                        ParallelParticleDriver(1, NTEILMAX + 1);
                    }
                    else
                    {
                        ParallelParticleDriver(1, NTEILMAX + 1);//20260408
                    }
                    if (Program.TransientReleaseEnabled) //20260521
                    {
                        int alive = 0;
                        for (int i = 1; i <= ActiveParticleMax; i++)
                        {
                            if (Program.ParticleSource[i] != 0 && Program.ParticleMass[i] > 0)
                            {
                                alive++;
                            }
                        }
                        Console.WriteLine($"[TRANSIENT_RELEASE] IWET={IWET} emissionFactor={releaseEmissionFactor:0.###} releasedThisStep={releasedThisStep} ActiveParticleMax={ActiveParticleMax} alive={alive}"); //20260521
                    }
                    if (Program.GasBuoyancyMode == 1)
                    {
                        PrintDensePhaseStats(IWET);
                    }
                    Console.WriteLine();

                    // Wait until conz4d file is written
                    if (ThreadWriteConz4dFile != null) // if Thread has been started -> wait until GFF--WriteThread has been finished
                    {
                        ThreadWriteConz4dFile.Join();  // wait, until thread has been finished
                        ThreadWriteConz4dFile = null;  // Release ressources				
                    }

                    //tranferring non-steady-state concentration fields
                    if (ISTATIONAER == Consts.TransientMode)
                    {
                        // Conz5d 是“当前时刻 transient 3D 场”，必须在 Conz5d->Conz4d 交换前写出 20260401
                        //double simTimeSeconds = Program.IWET * Program.TAUS;//20260402注释
                        double stepStartSimTimeSeconds = (Program.IWET - Program.SnapshotIWETBase) * Program.TAUS;//20260402
                        double simTimeSeconds = stepStartSimTimeSeconds + Program.TAUS;//20260402
                        bool isLastTransientStep = Program.IWET >= Program.MeteoTimeSer.Count;
                        TryWrite3DSnapshotFrame(WriteClass, simTimeSeconds, isLastTransientStep);
                        Console.WriteLine("[DBG] Enter TransferNonSteadyStateConcentrations");//20250106
                        TransferNonSteadyStateConcentrations(WriteClass, ref ThreadWriteConz4dFile);
                    }

                    if (FirstLoop)
                    {
                        ProgramWriters.LogfileGralCoreInfo(ResultFileZipped, GffFiles);
                        FirstLoop = false;
                    }

                    //Correction of concentration volume in each cell and for each gridded receptor 
                    VolumeCorrection();

                    //time needed for dispersion computations
                    double CalcTimeDispersion = (Environment.TickCount - StartTime) * 0.001 - CalcTimeWindfield;
                    Console.WriteLine();
                    Console.WriteLine("[PG] after DISPERSION Program.cs (common join)");
                    Console.WriteLine("Total simulation time [s]: " + (CalcTimeWindfield + CalcTimeDispersion).ToString("0.0"));
                    Console.WriteLine("Dispersion [s]: " + CalcTimeDispersion.ToString("0.0"));
                    Console.WriteLine("Flow field [s]: " + CalcTimeWindfield.ToString("0.0"));

                    if (LogLevel > Consts.LogLevelOff) // additional LOG-Output
                    {
                        LOG01_Output();
                    }

                    ZippedFile = IWET.ToString("00000") + ".grz";
                    try
                    {
                        if (File.Exists(ZippedFile))
                        {
                            File.Delete(ZippedFile); // delete existing files
                        }
                    }
                    catch { }

                    //output of 2-D concentration files (concentrations, deposition, odour-files)
                    Console.WriteLine(">>> LEAVE DISPERSION LOOP (transient) <<<");
                    Console.WriteLine("Writing result files in a background thread");
                    recentWeatherSituation = IWET;
                    ThreadWrite2DConcentrationFiles = new Thread(() => WriteClass.Write2DConcentrations(recentWeatherSituation, ZippedFile));
                    ThreadWrite2DConcentrationFiles.Start(); // start writing thread

                    //receptor concentrations
                    if (ISTATIONAER == Consts.TransientMode)
                    {
                        WriteClass.WriteReceptorTimeseries(0);
                    }
                    WriteClass.WriteReceptorConcentrations();

                    //microscale flow-field at receptors
                    WriteClass.WriteMicroscaleFlowfieldReceptors();

                } //skipped situation if no entry in meteopgt.all could be found in transient GRAL mode
            } // loop for all meteorological situations

            // Write summarized emission per source group and 3D Concentration file
            if (ISTATIONAER == Consts.TransientMode)
            {
                ProgramWriters WriteClass = new ProgramWriters();
                if (WriteVerticalConcentration)
                {
                    WriteClass.Write3DTextConcentrations();
                }
                Console.WriteLine();
                ProgramWriters.LogfileGralCoreWrite("");
                ProgramWriters.ShowEmissionRate();
            }

            if (ThreadWriteGffFiles != null) // if Thread has been started -> wait until GFF--WriteThread has been finished
            {
                ThreadWriteGffFiles.Join(); // wait, until thread has been finished
                ThreadWriteGffFiles = null; // Release ressources
            }

            if (ThreadWrite2DConcentrationFiles != null) // if Thread has been started -> wait until *.grz WriteThread has been finished
            {
                ThreadWrite2DConcentrationFiles.Join(); // wait, until thread has been finished
                ThreadWrite2DConcentrationFiles = null; // Release ressources
            }

            // Clean transient files
            Delete_Temp_Files();

            //Write receptor statistical error
            if (Program.ReceptorsAvailable)
            {
                ProgramWriters WriteClass = new ProgramWriters();
                WriteClass.WriteReceptorTimeseries(1);
            }

            ProgramWriters.LogfileGralCoreWrite("GRAL simulations finished at: " + DateTime.Now.ToString());
            ProgramWriters.LogfileGralCoreWrite(new String('-', 80));

            if (Program.IOUTPUT <= 0 && Program.WaitForConsoleKey) // not for Soundplan or no keystroke
            {
                Console.WriteLine();
                Console.WriteLine("GRAL simulations finished. Press any key to continue...");
                Program.CleanUpMemory();
                Console.ReadKey(true); 	// wait for a key input
            }

            Environment.Exit(0);        // Exit console
        }
        
        internal static long FlattenSnapshotVoxelKey(int i, int j, int k)
        {
            long ny = Program.NJJ + 2L;
            long nz = Program.NKK_Transient + 1L;
            return ((long)i * ny + j) * nz + k;
        }

        internal static void AccumulateIntraStepSnapshotFromConz5d(
            int reflexion_flag, float zcoord_nteil, float AHint, double masse, double Area_cart, float idtRaw,
            double xsi, double eta, int SG_nteil, double localStepEndSeconds)
        {
            if (reflexion_flag != 0) return;
            if (!Enable3DSnapshots || SnapshotSamplingMode != SnapshotSamplingModeKind.IntraStep) return;
            if (SnapshotIntervalSeconds <= 0 || idtRaw <= 0f) return;

            int indexK3d = TransientConcentration.BinarySearchTransient(zcoord_nteil - AHint);
            int indexI3d = (int)(xsi / Program.DXK) + 1;
            int indexJ3d = (int)(eta / Program.DYK) + 1;

            if (indexI3d < 1 || indexI3d > Program.NII + 1 ||
                indexJ3d < 1 || indexJ3d > Program.NJJ + 1 ||
                indexK3d < 1 || indexK3d > Program.NKK_Transient) return;

            double concContribution = masse * Program.GridVolume * Program.TAUS / (Area_cart * Program.DZK_Trans[indexK3d]);
            AccumulateIntraStepSnapshotCore(concContribution, localStepEndSeconds, idtRaw, indexI3d, indexJ3d, indexK3d);
        }

        internal static void AccumulateIntraStepSnapshotFromConz5dTransient(
            int reflexion_flag, float zcoord_nteil, float AHint, double mass_real, double Area_cart, float idtRaw,
            double xsi, double eta, int SG_nteil, double localStepEndSeconds)
        {
            if (reflexion_flag != 0) return;
            if (!Enable3DSnapshots || SnapshotSamplingMode != SnapshotSamplingModeKind.IntraStep) return;
            if (SnapshotIntervalSeconds <= 0 || idtRaw <= 0f) return;

            int indexK3d = TransientConcentration.BinarySearchTransient(zcoord_nteil - AHint);
            int indexI3d = (int)(xsi / Program.DXK) + 1;
            int indexJ3d = (int)(eta / Program.DYK) + 1;

            if (indexI3d < 1 || indexI3d > Program.NII + 1 ||
                indexJ3d < 1 || indexJ3d > Program.NJJ + 1 ||
                indexK3d < 1 || indexK3d > Program.NKK_Transient) return;

            double concContribution = mass_real * Program.TAUS / (Area_cart * Program.DZK_Trans[indexK3d]);
            AccumulateIntraStepSnapshotCore(concContribution, localStepEndSeconds, idtRaw, indexI3d, indexJ3d, indexK3d);
        }

        private static void AccumulateIntraStepSnapshotCore(
            double concContribution, double localStepEndSeconds, float idtRaw,
            int indexI3d, int indexJ3d, int indexK3d)
        {
            const double eps = 1e-12;
            if (Math.Abs(concContribution) <= eps) return;

            double stepStartTime = (Program.IWET - Program.SnapshotIWETBase) * Program.TAUS;
            double localStepStartSeconds = Math.Max(0.0, localStepEndSeconds - idtRaw);

            double globalStart = stepStartTime + localStepStartSeconds;
            double globalEnd = stepStartTime + localStepEndSeconds;

            if (globalEnd <= Program.SnapshotStartSeconds + eps) return;
            if (globalStart >= Program.SnapshotEndSeconds - eps) return;

            double start = Math.Max(globalStart, Program.SnapshotStartSeconds);
            double end = Math.Min(globalEnd, Program.SnapshotEndSeconds);
            if (end <= start + eps) return;

            double denom = Math.Max((double)idtRaw, eps); // fixed denominator = original idt
            long voxelKey = FlattenSnapshotVoxelKey(indexI3d, indexJ3d, indexK3d);

            long bucketId = (long)Math.Floor((start - Program.SnapshotStartSeconds) / Program.SnapshotIntervalSeconds);
            if (bucketId < 0) bucketId = 0;

            while (true)
            {
                double bucketStart = Program.SnapshotStartSeconds + bucketId * Program.SnapshotIntervalSeconds;
                if (bucketStart >= end - eps) break;

                double bucketEnd = bucketStart + Program.SnapshotIntervalSeconds;
                double overlap = Math.Min(end, bucketEnd) - Math.Max(start, bucketStart);

                if (overlap > eps)
                {
                    double part = concContribution * (overlap / denom);
                    var bucket = Program.SnapshotIntraBuckets.GetOrAdd(
                        bucketId,
                        _ => new System.Collections.Concurrent.ConcurrentDictionary<long, double>());
                    bucket.AddOrUpdate(voxelKey, part, (_, oldVal) => oldVal + part);
                }

                bucketId++;
            }
        }

        private static void TryWrite3DSnapshotFrame(ProgramWriters writeClass, double simTimeSeconds, bool forceFinalFrame)
        {
            if (ISTATIONAER != Consts.TransientMode || !Enable3DSnapshots) return;

            if (SnapshotSamplingMode == SnapshotSamplingModeKind.IntraStep)
            {
                TryWrite3DSnapshotFrameIntraStep(writeClass, simTimeSeconds, forceFinalFrame);
                return;
            }

            const double eps = 1e-9;
            if (simTimeSeconds + eps < SnapshotStartSeconds) return;

            bool withinEnd = simTimeSeconds <= SnapshotEndSeconds + eps;
            bool dueByInterval = withinEnd && (simTimeSeconds + eps >= SnapshotNextTimeSeconds);
            bool dueByFinal = forceFinalFrame && SnapshotWriteFinalFrame && withinEnd &&
                            (SnapshotLastWrittenTimeSeconds < 0 || Math.Abs(simTimeSeconds - SnapshotLastWrittenTimeSeconds) > eps);

            if (!dueByInterval && !dueByFinal) return;

            if (dueByInterval)
            {
                while (SnapshotNextTimeSeconds <= simTimeSeconds + eps)
                    SnapshotNextTimeSeconds += SnapshotIntervalSeconds;
            }

            int frameId = SnapshotFrameId++;
            if (writeClass.Write3DSnapshotConcentrationFrame(
                SnapshotOutputDir,
                frameId,
                simTimeSeconds,
                SnapshotUseSourceGroupSum,
                out string fileName))
            {
                SnapshotLastWrittenTimeSeconds = simTimeSeconds;
                ProgramWriters.Append3DSnapshotIndex(SnapshotOutputDir, frameId, simTimeSeconds, "StepEnd", fileName, false);
            }
        }

        private static void TryWrite3DSnapshotFrameIntraStep(ProgramWriters writeClass, double simTimeSeconds, bool forceFinalFrame)
        {
            const double eps = 1e-9;
            if (simTimeSeconds + eps < SnapshotStartSeconds) return;

            double closedHorizon = Math.Min(simTimeSeconds, SnapshotEndSeconds);
            double x = (closedHorizon - SnapshotStartSeconds) / SnapshotIntervalSeconds;
            long maxClosedBucketId = (long)Math.Floor(x - 1.0 + eps);

            while (SnapshotNextBucketId <= maxClosedBucketId)
            {
                long bucketId = SnapshotNextBucketId++;
                double frameTime = SnapshotStartSeconds + (bucketId + 1) * SnapshotIntervalSeconds;

                Program.SnapshotIntraBuckets.TryRemove(bucketId, out var bucketData);

                int frameId = SnapshotFrameId++;
                if (writeClass.Write3DSnapshotConcentrationFrameFromBucket(
                    SnapshotOutputDir,
                    frameId,
                    frameTime,
                    SnapshotUseSourceGroupSum,
                    bucketData,
                    out string fileName))
                {
                    SnapshotLastWrittenTimeSeconds = frameTime;
                    ProgramWriters.Append3DSnapshotIndex(SnapshotOutputDir, frameId, frameTime, "IntraStep", fileName, false);
                }
            }

            if (!forceFinalFrame || !SnapshotWriteFinalFrame) return;

            double finalTime = Math.Min(simTimeSeconds, SnapshotEndSeconds);
            if (finalTime + eps < SnapshotStartSeconds) return;

            bool alreadyWritten = SnapshotLastWrittenTimeSeconds >= 0 &&
                                Math.Abs(finalTime - SnapshotLastWrittenTimeSeconds) <= eps;

            double rel = (finalTime - SnapshotStartSeconds) / SnapshotIntervalSeconds;
            bool onBoundary = Math.Abs(rel - Math.Round(rel)) <= eps;
            if (alreadyWritten || onBoundary) return;

            long openBucketId = (long)Math.Floor(rel);
            if (openBucketId < SnapshotNextBucketId) return;

            Program.SnapshotIntraBuckets.TryRemove(openBucketId, out var openBucketData);

            int finalFrameId = SnapshotFrameId++;
            if (writeClass.Write3DSnapshotConcentrationFrameFromBucket(
                SnapshotOutputDir,
                finalFrameId,
                finalTime,
                SnapshotUseSourceGroupSum,
                openBucketData,
                out string finalName))
            {
                SnapshotLastWrittenTimeSeconds = finalTime;
                ProgramWriters.Append3DSnapshotIndex(SnapshotOutputDir, finalFrameId, finalTime, "IntraStep", finalName, true);
            }
        }


        private static void ReadDenseLightConfig(string cfgPath)
        {
            if (!File.Exists(cfgPath)) return;
            foreach (var raw in File.ReadAllLines(cfgPath))
            {
                var line = raw.Trim();
                if (line.Length == 0 || line.StartsWith("#")) continue;
                var eq = line.IndexOf('=');
                if (eq <= 0) continue;

                string key = line.Substring(0, eq).Trim();
                string val = line.Substring(eq + 1).Trim();
                int inlineComment = val.IndexOf('#'); //20260521 1852 allow timestamped inline comments in densegas.cfg values
                if (inlineComment >= 0) val = val.Substring(0, inlineComment).Trim(); //20260521 1852 parse value before inline comment

                // bools
                if (key.Equals("UseDenseGasMicroSpray", StringComparison.OrdinalIgnoreCase))
                    Program.UseDenseGasMicroSpray = val.Equals("true", StringComparison.OrdinalIgnoreCase);
                else if (key.Equals("UseDenseHorzRandomWalk", StringComparison.OrdinalIgnoreCase))
                    Program.UseDenseHorzRandomWalk = val.Equals("true", StringComparison.OrdinalIgnoreCase);
                else if (key.Equals("UseExplicitExitVelocityDirection", StringComparison.OrdinalIgnoreCase)) //20260521 1852 optional Trial 8 point-source velocity-vector switch
                    Program.UseExplicitExitVelocityDirection = val.Equals("true", StringComparison.OrdinalIgnoreCase) || val.Equals("1", StringComparison.OrdinalIgnoreCase);

                // ints
                else if (key.Equals("GasBuoyancyMode", StringComparison.OrdinalIgnoreCase))
                    int.TryParse(val, out Program.GasBuoyancyMode);

                // floats (Culture-invariant)
                else if (key.Equals("GasMolWeight", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.GasMolWeight);
                else if (key.Equals("GasTempK", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.GasTempK);
                else if (key.Equals("AmbientTempK", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.AmbientTempK);
                else if (key.Equals("GasDensityRatio0", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.GasDensityRatio0);
                else if (key.Equals("InitialCloudDepth", StringComparison.OrdinalIgnoreCase))
                {
                    if (float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out var d))
                    {
                        _initCloudDepth = MathF.Max(0.1f, d);
                    }
                }
                else if (key.Equals("InitialCloudRadius", StringComparison.OrdinalIgnoreCase))
                {
                    if (float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out var r))
                    {
                        _initCloudRadius = MathF.Max(0.1f, r);
                    }
                }
                else if (key.Equals("DenseVertCoeff", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.DenseVertCoeff);
                //20260304
                else if (key.Equals("DenseHorzCoeff", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.DenseHorzCoeff);
                else if (key.Equals("DenseHorzRandomMaxStepFrac", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.DenseHorzRandomMaxStepFrac);
                else if (key.Equals("DenseEpsNeutral", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.DenseEpsNeutral);
                else if (key.Equals("DenseMixFactor", StringComparison.OrdinalIgnoreCase))
                {
                    if (float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out var mix))
                    {
                        Program.DenseMixFactor = MathF.Max(0.0f, mix);
                    }
                }
                // 20260510 dense-stage spatial entrainment multiplier
                else if (key.Equals("DenseEntrainmentFactor", StringComparison.OrdinalIgnoreCase))
                {
                    if (float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out var entrainment))
                    {
                        Program.DenseEntrainmentFactor = MathF.Max(0.0f, entrainment);
                    }
                }
                else if (key.Equals("DenseGroundEps", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.DenseGroundEps);
                else if (key.Equals("DenseSpreadFactor", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.DenseSpreadFactor);
                else if (key.Equals("DenseUgMax", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.DenseUgMax);
                else if (key.Equals("DenseReflA", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.DenseReflA);
                else if (key.Equals("DenseReflB", StringComparison.OrdinalIgnoreCase))
                    float.TryParse(val, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out Program.DenseReflB);
                else if (key.Equals("DenseColumnSmooth3x3", StringComparison.OrdinalIgnoreCase))
                    Program.DenseColumnSmooth3x3 = val.Equals("true", StringComparison.OrdinalIgnoreCase);
            }
        }
        private static void InitializeDenseColumnFields()
        {
            int nx = Math.Max(2, NII + 2);
            int ny = Math.Max(2, NJJ + 2);

            DenseColumnSumH = new double[nx][];
            DenseColumnSumRho = new double[nx][];
            DenseColumnCount = new int[nx][];
            DenseColumnHbulk = new float[nx][];
            DenseColumnRhoBulk = new float[nx][];
            DenseColumnBbulk = new float[nx][];

            for (int i = 0; i < nx; i++)
            {
                DenseColumnSumH[i] = new double[ny];
                DenseColumnSumRho[i] = new double[ny];
                DenseColumnCount[i] = new int[ny];
                DenseColumnHbulk[i] = new float[ny];
                DenseColumnRhoBulk[i] = new float[ny];
                DenseColumnBbulk[i] = new float[ny];
            }

            DenseColumnLocks = new object[256];
            for (int i = 0; i < DenseColumnLocks.Length; i++)
            {
                DenseColumnLocks[i] = new object();
            }
        }

        private static void ResetDenseColumnFields()
        {
            if (DenseColumnCount == null || DenseColumnCount.Length == 0)
            {
                return;
            }

            for (int i = 0; i < DenseColumnCount.Length; i++)
            {
                Array.Clear(DenseColumnSumH[i], 0, DenseColumnSumH[i].Length);
                Array.Clear(DenseColumnSumRho[i], 0, DenseColumnSumRho[i].Length);
                Array.Clear(DenseColumnCount[i], 0, DenseColumnCount[i].Length);
                Array.Clear(DenseColumnHbulk[i], 0, DenseColumnHbulk[i].Length);
                Array.Clear(DenseColumnRhoBulk[i], 0, DenseColumnRhoBulk[i].Length);
                Array.Clear(DenseColumnBbulk[i], 0, DenseColumnBbulk[i].Length);
            }
        }
        
        //新增N/D统计函数 20260312
        private static void PrintDensePhaseStats(int iwet)
        {
            if (GasBuoyancyMode != 1 || DenseStates == null || DenseStates.Length <= 1) return;

            int imax = Math.Min(NTEILMAX, DenseStates.Length - 1);
            int nNeutral = 0;
            int nDenseActive = 0;
            int nAlive = 0;

            float rhoAir = 1.2041f * (273.15f / MathF.Max(200f, AmbientTempK));

            double sumRhoBulk = 0.0;   // dense-active 粒子的平均 bulk 密度
            double sumDeltaAbs = 0.0;  // (rhoBulk - rhoAir) 的平均绝对差
            double sumDeltaRel = 0.0;  // (rhoBulk/rhoAir - 1) 的平均相对差

            for (int ip = 1; ip <= imax; ip++)
            {
                // active particles in current run
                if (ParticleSource[ip] == 0 || ParticleMass[ip] <= 0) continue;

                nAlive++;
                if (DenseStates[ip].phase == DensePhase.Neutral) nNeutral++;
                else nDenseActive++;
            
            // 优先用局地柱平均（若有），否则回退到粒子态 rhoBulk
            float rho = (DenseStates[ip].RhoBulkLocal > 0f) ? DenseStates[ip].RhoBulkLocal : DenseStates[ip].rhoBulk;
            if (rho <= 0f) rho = rhoAir;

            float dAbs = rho - rhoAir;
            float dRel = dAbs / MathF.Max(1e-6f, rhoAir);

            sumRhoBulk += rho;
            sumDeltaAbs += dAbs;
            sumDeltaRel += dRel;
            }
            float ratioNeutralTotal = (float)nNeutral / Math.Max(1, NTEILMAX);
            float ratioNeutralAlive = (float)nNeutral / Math.Max(1, nAlive);
            float avgRhoBulk = (nDenseActive > 0) ? (float)(sumRhoBulk / nDenseActive) : 0f;
            float avgDeltaAbs = (nDenseActive > 0) ? (float)(sumDeltaAbs / nDenseActive) : 0f;
            float avgDeltaRel = (nDenseActive > 0) ? (float)(sumDeltaRel / nDenseActive) : 0f;
            float avgRelDensity = (rhoAir > 0f) ? avgRhoBulk / rhoAir : 0f;

            Console.WriteLine(
                $"[DENSE_STAT] IWET={iwet} " +
                $"Nneutral={nNeutral} Ndense-active={nDenseActive} Alive={nAlive} " +
                $"Neutral/Total={ratioNeutralTotal:P2} Neutral/Alive={ratioNeutralAlive:P2}" +
                $"rhoAir={rhoAir:F4} rhoBulkAvg={avgRhoBulk:F4} " +
                $"relRhoAvg={avgRelDensity:F4} dRhoAbsAvg={avgDeltaAbs:F4} dRhoRelAvg={avgDeltaRel:E3}");
        }

        private static float _initCloudDepth = 1.0f;
        private static float _initCloudRadius = 1.0f;

        public static bool DebugDense = true;
        public static int DebugDensePrintCount = 0;
    }
}



