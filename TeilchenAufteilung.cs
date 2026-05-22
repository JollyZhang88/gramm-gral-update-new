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

namespace GRAL_2001
{
    class ParticleManagement
    {
        /// <summary>
        /// Calculate the Sum of all Emission Rates 
        /// Assign the particles to the sources
        /// </summary>
        public static int Calculate()
        {
            //total emissions (all sources)
            Program.PS_PartNumb.Initialize(); Program.PS_PartSum = 0;
            Program.LS_PartNumb.Initialize(); Program.LS_PartSum = 0;
            Program.AS_PartNumb.Initialize(); Program.AS_PartSum = 0;
            Program.TS_PartNumb.Initialize(); Program.TS_PartSum = 0;
            Program.VS_PartNumb.Initialize(); Program.VS_PartSum = 0;//20260310 add volume sources

            double sum_emission = 0.0; // Sum for all sources
            int Sum_of_Particles = 0; // sum for alternative approach

            // Calculate sum of all emissions
            for (int i = 1; i <= Program.PS_Count; i++) // Point sources
            {
                sum_emission += Program.PS_ER[i];
            }

            for (int i = 1; i <= Program.LS_Count; i++) // Line sources
            {
                sum_emission += Program.LS_ER[i];
            }

            for (int i = 1; i <= Program.TS_Count; i++) // Portals
            {
                sum_emission += Program.TS_ER[i];
            }

            for (int i = 1; i <= Program.AS_Count; i++) // Area sources
            {
                sum_emission += Program.AS_ER[i];
            }

            for (int i = 1; i <= Program.VS_Count; i++) // Volume sources
            {
                sum_emission += Program.VS_ER[i];
            }

            int PS_Min_Particles = Set_Min_Particles_PS_TS(Program.PS_Count);
            int LS_Min_Particles = Set_Min_Particles_LS(Program.LS_Count);
            int TS_Min_Particles = Set_Min_Particles_PS_TS(Program.TS_Count);
            int AS_Min_Particles = Set_Min_Particles_AS(Program.AS_Count);
            int VS_Min_Particles = Set_Min_Particles_AS(Program.VS_Count);

            // --- DEBUG & SAFETY GUARD: sum_emission and NTEILMAX ---
            if (double.IsNaN(sum_emission) || double.IsInfinity(sum_emission) || sum_emission <= 0.0)
            {
                Console.WriteLine($"[ERR] sum_emission invalid: {sum_emission}. " +
                                  $"PS_Count={Program.PS_Count}, LS_Count={Program.LS_Count}, TS_Count={Program.TS_Count}, AS_Count={Program.AS_Count}");
                return -1; // 直接返回，让上层知道错误
            }
            if (Program.NTEILMAX <= 0)
            {
                Console.WriteLine($"[ERR] Program.NTEILMAX invalid: {Program.NTEILMAX}");
                return -1;
            }
            Console.WriteLine($"[OK] Emission sum={sum_emission:0.###}, NTEILMAX={Program.NTEILMAX} -> start particle assignment...");
            Console.WriteLine("[OK] particle assignment finished.");

            // 原来是：unit = (float)Program.NTEILMAX / sum_emission;
            double unit = Program.NTEILMAX / sum_emission;

            /*
            //assignment of particles for each source
            double unit = (float)Program.NTEILMAX / sum_emission;
            */

            int sum = 0;
            for (int i = 1; i <= Program.PS_Count; i++)
            {
                Program.PS_PartNumb[i] = Convert.ToInt32(unit * Program.PS_ER[i]);

                if (Program.PS_PartNumb[i] < PS_Min_Particles)
                {
                    Program.PS_PartNumb[i] = PS_Min_Particles;
                }

                Program.PS_PartSum += Program.PS_PartNumb[i];
                sum += Program.PS_PartNumb[i];
                Sum_of_Particles += Program.PS_PartNumb[i];

                if (Program.LogLevel == Consts.LogLevelRefPart) // Show all particle-numbers in the case of Log-Level 02
                {
                    Console.WriteLine("PS " + i.ToString() + " : " + Math.Abs(Program.PS_PartNumb[i]).ToString());
                }
            }

            for (int i = 1; i <= Program.LS_Count; i++)
            {
                Program.LS_PartNumb[i] = Convert.ToInt32(unit * Program.LS_ER[i]);

                if (Program.LS_PartNumb[i] < LS_Min_Particles) // avoid Null particle sources
                {
                    Program.LS_PartNumb[i] = LS_Min_Particles;
                }

                Program.LS_PartSum += Program.LS_PartNumb[i];
                sum += Program.LS_PartNumb[i];
                Sum_of_Particles += Program.LS_PartNumb[i];

                if (Program.LogLevel == Consts.LogLevelRefPart) // Show all particle-numbers in the case of Log-Level 02
                {
                    Console.WriteLine("LS " + i.ToString() + " : " + Math.Abs(Program.LS_PartNumb[i]).ToString());
                }
            }

            for (int i = 1; i <= Program.TS_Count; i++)
            {
                Program.TS_PartNumb[i] = Convert.ToInt32(unit * Program.TS_ER[i]);

                if (Program.TS_PartNumb[i] < TS_Min_Particles)
                {
                    Program.TS_PartNumb[i] = TS_Min_Particles;
                }

                Program.TS_PartSum += Program.TS_PartNumb[i];
                sum += Program.TS_PartNumb[i];
                Sum_of_Particles += Program.TS_PartNumb[i];

                if (Program.LogLevel == Consts.LogLevelRefPart) // Show all particle-numbers in the case of Log-Level 02
                {
                    Console.WriteLine("TS " + i.ToString() + " : " + Math.Abs(Program.TS_PartNumb[i]).ToString());
                }
            }

            for (int i = 1; i <= Program.AS_Count; i++)
            {
                Program.AS_PartNumb[i] = Convert.ToInt32(unit * Program.AS_ER[i]);

                if (Program.AS_PartNumb[i] < AS_Min_Particles) // avoid Null particle sources
                {
                    Program.AS_PartNumb[i] = AS_Min_Particles;
                }

                Program.AS_PartSum += Program.AS_PartNumb[i];
                sum += Program.AS_PartNumb[i];
                Sum_of_Particles += Program.AS_PartNumb[i];

                if (Program.LogLevel == Consts.LogLevelRefPart) // Show all particle-numbers in the case of Log-Level 02
                {
                    Console.WriteLine("AS " + i.ToString() + " : " + Math.Abs(Program.AS_PartNumb[i]).ToString());
                }
            }
            
            for (int i = 1; i <= Program.VS_Count; i++)
            {
                double emWeight = GetVSEmissionWeight(i);
                Program.VS_PartNumb[i] = Convert.ToInt32(unit * emWeight);

                if (Program.VS_PartNumb[i] < VS_Min_Particles)
                {
                    Program.VS_PartNumb[i] = VS_Min_Particles;
                }

                Program.VS_PartSum += Program.VS_PartNumb[i];
                sum += Program.VS_PartNumb[i];
                Sum_of_Particles += Program.VS_PartNumb[i];

                if (Program.LogLevel == Consts.LogLevelRefPart)
                {
                    Console.WriteLine("VS " + i.ToString() + " : " + Math.Abs(Program.VS_PartNumb[i]).ToString());
                }
            }


            Program.ParticleMassMean = sum_emission / 3600 / Program.TPS * 1000000000 / Program.GridVolume / Program.TAUS;

            string err = "Using a total of " + Sum_of_Particles.ToString() + " particles (" + Math.Round((Sum_of_Particles / Program.TAUS), 0).ToString() + " part./s) within the model domain";
            Console.WriteLine("[PG] HIT A — after 'Using a total of ...'");
            Console.WriteLine(err);
            ProgramWriters.LogfileGralCoreWrite(err);
            ProgramWriters.LogfileGralCoreWrite("");

            Program.ParticlesPerReleaseStep = Sum_of_Particles; //20260521
            int capacity = Sum_of_Particles; //20260521
            if (Program.TransientReleaseEnabled) //20260521
            {
                int weatherSteps = Math.Max(1, Program.MeteoTimeSer.Count);
                long requestedCapacity = (long)Sum_of_Particles * weatherSteps;
                if (requestedCapacity > int.MaxValue - 2)
                {
                    Console.WriteLine($"[ERR] transient particle capacity too large: perStep={Sum_of_Particles}, weatherSteps={weatherSteps}");
                    return -1;
                }

                capacity = (int)requestedCapacity;
                Program.TotalParticleCapacity = capacity;
                Program.ActiveParticleMax = 0;

                string releaseInfo = $"Transient continuous release: particles per release step={Program.ParticlesPerReleaseStep}, weather steps={weatherSteps}, total particle capacity={Program.TotalParticleCapacity}";
                Console.WriteLine(releaseInfo);
                ProgramWriters.LogfileGralCoreWrite(releaseInfo);
            }
            else
            {
                Program.TotalParticleCapacity = capacity; //20260521
                Program.ActiveParticleMax = 0; //20260521
            }

            // Resize Arrays 
            Program.ParticleSource = GC.AllocateUninitializedArray<Int32>(capacity + 1); //20260521
            Program.ParticleSG = GC.AllocateUninitializedArray<byte>(capacity + 1); //20260521
            Program.Xcoord = GC.AllocateUninitializedArray<double>(capacity + 1); //20260521
            Program.YCoord = GC.AllocateUninitializedArray<double>(capacity + 1); //20260521
            Program.ZCoord = GC.AllocateUninitializedArray<float>(capacity + 1); //20260521
            Program.ParticleMass = GC.AllocateUninitializedArray<double>(capacity + 1); //20260521
            Program.SourceType = GC.AllocateUninitializedArray<byte>(capacity + 1); //20260521
            Program.ParticleVsed = GC.AllocateUninitializedArray<float>(capacity + 1);         // sedimentation velocity of one lagrangian particle //20260521
            Program.ParticleVdep = GC.AllocateUninitializedArray<float>(capacity + 1);         // deposition velocity of one lagrangian particle //20260521
            Program.ParticleMode = GC.AllocateUninitializedArray<byte>(capacity + 1); //20260521

            Console.WriteLine();

            return capacity; //20260521
        }

        private static int Set_Min_Particles_PS_TS(int source_count)
        {
            int Min_Particles = 10;
            if (source_count < 2000)
            {
                Min_Particles = 20;
            }
            else if (source_count > 30000)
            {
                Min_Particles = 5;
            }

            return Min_Particles;
        }

        private static int Set_Min_Particles_LS(int source_count)
        {
            int Min_Particles = 5;
            if (source_count < 2000)
            {
                Min_Particles = 16;
            }
            else if (source_count < 10000)
            {
                Min_Particles = 8;
            }
            return Min_Particles;
        }

        private static int Set_Min_Particles_AS(int source_count)
        {
            int Min_Particles = 2;
            if (source_count < 2000)
            {
                Min_Particles = 6;
            }
            else if (source_count < 10000)
            {
                Min_Particles = 4;
            }

            return Min_Particles;
        }
        private static double GetVSEmissionWeight(int i)
        {
            if (Program.VS_Instant != null && i < Program.VS_Instant.Length && Program.VS_Instant[i])
            {
                double totalKg = (Program.VS_TotalRelease != null && i < Program.VS_TotalRelease.Length)
                    ? Math.Max(0.0, Program.VS_TotalRelease[i])
                    : 0.0;

                // 折算成 kg/h，用于与其它源统一分配权重
                return totalKg * 3600.0 / Math.Max(1e-6, Program.TAUS);
            }

            return (Program.VS_ER != null && i < Program.VS_ER.Length) ? Math.Max(0.0, Program.VS_ER[i]) : 0.0;
        }

    }
}
