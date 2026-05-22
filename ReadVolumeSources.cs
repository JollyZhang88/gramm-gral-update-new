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
using System.Collections.Generic;
using System.IO;

namespace GRAL_2001
{
    class ReadVolumeSources
    {
        /// <summary>
        /// Read Volume Sources from the file "volume.dat"
        /// </summary>
    	public static void Read()
        {

            List<SourceData> VQ = new List<SourceData>();
            VQ.Add(new SourceData()); // index 0 dummy

            double totalemission = 0.0;
            int countrealsources = 0;
            double[] emission_sourcegroup = new double[101];

            int lineNo = 1;
            using StreamReader read = new StreamReader("volume.dat");

            Deposition Dep = new Deposition();

            try
            {
                string line = read.ReadLine(); // header
                while ((line = read.ReadLine()) != null)
                {
                    lineNo++;
                    if (string.IsNullOrWhiteSpace(line)) continue;

                    string[] text = line.Split(new char[] { ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                    if (text.Length < 7) continue; // minimal columns

                    double x = Convert.ToDouble(text[0].Replace(".", Program.Decsep));
                    double y = Convert.ToDouble(text[1].Replace(".", Program.Decsep));
                    float z = Convert.ToSingle(text[2].Replace(".", Program.Decsep));
                    double r = Convert.ToDouble(text[3].Replace(".", Program.Decsep));
                    float h = Convert.ToSingle(text[4].Replace(".", Program.Decsep));
                    double er = Convert.ToDouble(text[5].Replace(".", Program.Decsep));
                    short sg = Convert.ToInt16(text[6]);
                    
                    if (r <= 0 || h <= 0) continue;
                    // XY domain check using cylinder footprint
                    double xsiMin = (x - r) - Program.IKOOAGRAL;
                    double xsiMax = (x + r) - Program.IKOOAGRAL;
                    double etaMin = (y - r) - Program.JKOOAGRAL;
                    double etaMax = (y + r) - Program.JKOOAGRAL;
                    if (!(etaMin > Program.EtaMinGral && xsiMin > Program.XsiMinGral &&
                        etaMax < Program.EtaMaxGral && xsiMax < Program.XsiMaxGral))
                    {
                        continue;
                    }
                    
                    int sgIndex = Program.Get_Internal_SG_Number(sg);
                    if (sgIndex < 0) continue;

                    SourceData sd = new SourceData
                    {
                        X1 = x,
                        Y1 = y,
                        Z1 = z,
                        X2 = r,     // store radius
                        Z2 = h,     // store height
                        ER = er,
                        SG = sg,
                        Mode = 0,
                        Vdep = 0f,
                        Vsed = 0f,
                        ER_dep = 0.0
                    };
                    
                    // optional deposition columns
                    if (text.Length > 7)
                    {
                        sd.Mode = Convert.ToByte(text[7]);
                    }
                    if (text.Length > 8)
                    {
                        sd.Vdep = Convert.ToSingle(text[8].Replace(".", Program.Decsep));
                    }
                    if (text.Length > 9)
                    {
                        sd.Vsed = Convert.ToSingle(text[9].Replace(".", Program.Decsep));
                    }
                    if (text.Length > 10)
                    {
                        sd.ER_dep = Convert.ToDouble(text[10].Replace(".", Program.Decsep));
                    }

                    VQ.Add(sd);
                    totalemission += sd.ER;
                    emission_sourcegroup[sgIndex] += sd.ER;
                    countrealsources++;
                }
            }
            catch
            {
                string err = "Error when reading file volume.dat in line " +  lineNo.ToString();
                Console.WriteLine(err + " Execution stopped: press ESC to stop");
                ProgramWriters.LogfileProblemreportWrite(err);

                if (Program.IOUTPUT <= 0 && Program.WaitForConsoleKey) // not for Soundplan or no keystroke
                {
                    Program.CleanUpMemory();
                    while (!(Console.KeyAvailable && Console.ReadKey(true).Key == ConsoleKey.Escape))
                    {
                        ;
                    }
                }

                Environment.Exit(0);
            }

            int counter = VQ.Count + 1;
            Program.VS_Count = VQ.Count - 1;

            // allocate arrays
            Program.VS_ER = GC.AllocateUninitializedArray<double>(counter);
            Program.VS_X = GC.AllocateUninitializedArray<double>(counter);
            Program.VS_Y = GC.AllocateUninitializedArray<double>(counter);
            Program.VS_Z = GC.AllocateUninitializedArray<float>(counter);
            Program.VS_R = GC.AllocateUninitializedArray<float>(counter);
            Program.VS_H = GC.AllocateUninitializedArray<float>(counter);

            Program.VS_SG = GC.AllocateUninitializedArray<byte>(counter);
            Program.VS_PartNumb = GC.AllocateUninitializedArray<int>(counter);
            Program.VS_Mode = GC.AllocateUninitializedArray<byte>(counter);
            Program.VS_V_Dep = GC.AllocateUninitializedArray<float>(counter);
            Program.VS_V_sed = GC.AllocateUninitializedArray<float>(counter);
            Program.VS_ER_Dep = GC.AllocateUninitializedArray<float>(counter);
            Program.VS_Absolute_Height = GC.AllocateUninitializedArray<bool>(counter);

            // if these exist in ProgramDeclarations, initialize them too
            Program.VS_Shape = GC.AllocateUninitializedArray<byte>(counter);
            Program.VS_ZRef = GC.AllocateUninitializedArray<byte>(counter);
            Program.VS_GroundAttach = GC.AllocateUninitializedArray<bool>(counter);
            Program.VS_Instant = GC.AllocateUninitializedArray<bool>(counter);
            Program.VS_ReleaseTimeS = GC.AllocateUninitializedArray<float>(counter);
            Program.VS_Released = GC.AllocateUninitializedArray<bool>(counter);

            // optional compatibility fields
            Program.VS_DX = GC.AllocateUninitializedArray<float>(counter);
            Program.VS_DY = GC.AllocateUninitializedArray<float>(counter);
            Program.VS_DZ = GC.AllocateUninitializedArray<float>(counter);

            for (int i = 1; i < VQ.Count; i++)
            {
                Program.VS_X[i] = VQ[i].X1;
                Program.VS_Y[i] = VQ[i].Y1;

                Program.VS_Absolute_Height[i] = (VQ[i].Z1 < 0);
                Program.VS_Z[i] = Math.Abs(VQ[i].Z1);

                Program.VS_R[i] = (float)VQ[i].X2;
                Program.VS_H[i] = VQ[i].Z2;

                Program.VS_DX[i] = 2f * Program.VS_R[i];
                Program.VS_DY[i] = 2f * Program.VS_R[i];
                Program.VS_DZ[i] = Program.VS_H[i];

                Program.VS_ER[i] = VQ[i].ER;
                Program.VS_SG[i] = (byte)VQ[i].SG;
                Program.VS_Mode[i] = VQ[i].Mode;
                Program.VS_V_Dep[i] = VQ[i].Vdep;
                Program.VS_V_sed[i] = VQ[i].Vsed;
                Program.VS_ER_Dep[i] = (float)VQ[i].ER_dep;

                Program.VS_Shape[i] = 2;          // cylinder
                Program.VS_ZRef[i] = 0;           // center
                Program.VS_GroundAttach[i] = false;
                Program.VS_Instant[i] = false;
                Program.VS_ReleaseTimeS[i] = 0f;
                Program.VS_Released[i] = false;
            }

            string info = "Total number of volume sources: " + countrealsources.ToString();
            Console.WriteLine(info);
            ProgramWriters.LogfileGralCoreWrite(info);

            string unit = Program.Odour ? "[MOU/h]: " : "[kg/h]: ";
            info = "Total volume emission " + unit + totalemission.ToString("0.000");
            Console.Write(info);
            ProgramWriters.LogfileGralCoreWrite(info);

            Console.Write(" (");
            for (int im = 0; im < Program.SourceGroups.Count; im++)
            {
                info = "  SG " + Program.SourceGroups[im] + unit + emission_sourcegroup[im].ToString("0.000");
                Console.Write(info);
                ProgramWriters.LogfileGralCoreWrite(info);
            }
            Console.WriteLine(" )");

            VQ.Clear();
            VQ.TrimExcess();
        }
    }
}
