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
    /// <summary>
    /// Read Point Sources
    /// </summary>
    class ReadPointSources
    {
        /// <summary>
        /// Read the point source data from the file "point.dat" and create and fill the point source arrays
        /// </summary>
        public static void Read()
        {

            List<SourceData> PQ = new List<SourceData>();

            double totalemission = 0;
            int countrealsources = 0;
            double[] emission_sourcegroup = new double[101];

            PQ.Add(new SourceData());

            Deposition Dep = new Deposition();

            StreamReader read = new StreamReader("point.dat");
            try
            {
                string[] text = new string[1];
                string text1;
                //text1 = read.ReadLine();//20260408 注释
                //text1 = read.ReadLine();//20260408 注释
                //20260408 读取标题行和表头行，判断是否包含点源速度方向列
                _ = read.ReadLine(); // title line
                string headerLine = read.ReadLine() ?? string.Empty;

                bool hasDirCols =
                    headerLine.IndexOf("DirX", StringComparison.OrdinalIgnoreCase) >= 0 &&
                    headerLine.IndexOf("DirY", StringComparison.OrdinalIgnoreCase) >= 0 &&
                    headerLine.IndexOf("DirZ", StringComparison.OrdinalIgnoreCase) >= 0;

                while ((text1 = read.ReadLine()) != null)
                {
                    text = text1.Split(new char[] { ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                    //if (text.Length > 9)
                    if (text.Length > 10)
                    {
                        double xsi = Convert.ToDouble(text[0].Replace(".", Program.Decsep)) - Program.IKOOAGRAL;
                        double eta = Convert.ToDouble(text[1].Replace(".", Program.Decsep)) - Program.JKOOAGRAL;
                        double dia = Convert.ToDouble(text[8].Replace(".", Program.Decsep));

                        //excluding all point sources outside GRAL domain
                        if (((eta - dia * 0.5) > Program.EtaMinGral) && ((xsi - dia * 0.5) > Program.XsiMinGral) && ((eta + dia * 0.5) < Program.EtaMaxGral) && ((xsi + dia * 0.5) < Program.XsiMaxGral))
                        {
                            //excluding all point sources with undesired source groups
                            {
                                Int16 SG = Convert.ToInt16(text[10]);
                                int SG_index = Program.Get_Internal_SG_Number(SG); // get internal SG number

                                if (SG_index >= 0)
                                {
                                    SourceData sd = new SourceData();
                                    sd.X1 = Convert.ToDouble(text[0].Replace(".", Program.Decsep));
                                    sd.Y1 = Convert.ToDouble(text[1].Replace(".", Program.Decsep));
                                    sd.Z1 = Convert.ToSingle(text[2].Replace(".", Program.Decsep));
                                    
                                    if (countrealsources < 5)//20260202
                                    {
                                        Console.WriteLine($"[DBG] PS read raw z={sd.Z1} x={sd.X1} y={sd.Y1}");
                                    }

                                    sd.ER = Convert.ToDouble(text[3].Replace(".", Program.Decsep));
                                    sd.V = Convert.ToSingle(text[7].Replace(".", Program.Decsep));
                                    sd.D = Convert.ToSingle(text[8].Replace(".", Program.Decsep));
                                    sd.T = Convert.ToSingle(text[9].Replace(".", Program.Decsep));
                                    sd.SG = Convert.ToInt16(text[10]);
                                    //20260408 point source velocity direction x/y/z
                                    int depStartIndex = 11; // old format default
                                    sd.VDirX = 0f;
                                    sd.VDirY = 0f;
                                    sd.VDirZ = 1f;

                                    if (hasDirCols)
                                    {
                                        if (text.Length <= 13) throw new IOException();

                                        sd.VDirX = Convert.ToSingle(text[11].Replace(".", Program.Decsep));
                                        sd.VDirY = Convert.ToSingle(text[12].Replace(".", Program.Decsep));
                                        sd.VDirZ = Convert.ToSingle(text[13].Replace(".", Program.Decsep));

                                        float n = MathF.Sqrt(sd.VDirX * sd.VDirX + sd.VDirY * sd.VDirY + sd.VDirZ * sd.VDirZ);
                                        if (n > 1e-6f)
                                        {
                                            sd.VDirX /= n;
                                            sd.VDirY /= n;
                                            sd.VDirZ /= n;
                                        }
                                        else
                                        {
                                            sd.VDirX = 0f;
                                            sd.VDirY = 0f;
                                            sd.VDirZ = 1f;
                                        }

                                        depStartIndex = 14;
                                    }
                                    sd.Mode = 0; // standard mode = concentration only
                                    totalemission += sd.ER;
                                    emission_sourcegroup[SG_index] += sd.ER;
                                    countrealsources++;

                                    sd.TimeSeriesTemperature = GetTransientTimeSeriesIndex.GetIndex(Program.PS_TimeSerTempValues, "Temp@_", text);
                                    sd.TimeSeriesVelocity = GetTransientTimeSeriesIndex.GetIndex(Program.PS_TimeSerVelValues, "Vel@_", text);

                                    //if (text.Length > 15) // deposition data available
                                    if (text.Length >= depStartIndex + 7)
                                    {
                                        //Dep.Dep_Start_Index = 11; // start index for point sources
                                        Dep.Dep_Start_Index = depStartIndex; // start index for point sources   20260408
                                        Dep.SD = sd;
                                        Dep.SourceData = PQ;
                                        Dep.Text = text;
                                        if (Dep.Compute() == false)
                                        {
                                            throw new IOException();
                                        }
                                    }
                                    else // no depositon
                                    {
                                        PQ.Add(sd);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            catch
            {
                string err = "Error when reading file point.dat in line " + (countrealsources + 3).ToString() + " Execution stopped: press ESC to stop";
                Console.WriteLine(err);
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
            read.Close();
            read.Dispose();

            int counter = PQ.Count + 1;
            Program.PS_Count += PQ.Count - 1;

            // Copy Lists to global arrays
            Program.PS_ER = GC.AllocateUninitializedArray<double>(counter);
            Program.PS_X = GC.AllocateUninitializedArray<double>(counter);
            Program.PS_Y = GC.AllocateUninitializedArray<double>(counter);
            Program.PS_Z = GC.AllocateUninitializedArray<float>(counter);
            Program.PS_V = GC.AllocateUninitializedArray<float>(counter);
            Program.PS_T = GC.AllocateUninitializedArray<float>(counter);
            Program.PS_D = GC.AllocateUninitializedArray<float>(counter);
            Program.PS_SG = GC.AllocateUninitializedArray<byte>(counter);
            Program.PS_PartNumb = GC.AllocateUninitializedArray<int>(counter);
            Program.PS_Mode = GC.AllocateUninitializedArray<byte>(counter);
            Program.PS_V_Dep = GC.AllocateUninitializedArray<float>(counter);
            Program.PS_V_sed = GC.AllocateUninitializedArray<float>(counter);
            Program.PS_ER_Dep = GC.AllocateUninitializedArray<float>(counter);
            Program.PS_Absolute_Height = GC.AllocateUninitializedArray<bool>(counter);
            Program.PS_TimeSeriesTemperature = GC.AllocateUninitializedArray<int>(counter);
            Program.PS_TimeSeriesVelocity = GC.AllocateUninitializedArray<int>(counter);
            Program.PS_VDirX = GC.AllocateUninitializedArray<float>(counter);//20260408 point source velocity direction x
            Program.PS_VDirY = GC.AllocateUninitializedArray<float>(counter);//20260408 point source velocity direction y
            Program.PS_VDirZ = GC.AllocateUninitializedArray<float>(counter);//20260408 point source velocity direction z


            for (int i = 1; i < PQ.Count; i++)
            {
                Program.PS_ER[i] = PQ[i].ER;
                Program.PS_X[i] = PQ[i].X1;
                Program.PS_Y[i] = PQ[i].Y1;

                if (PQ[i].Z1 < 0) // negative value = absolute height
                {
                    Program.PS_Absolute_Height[i] = true;
                    if (Program.Topo != Consts.TerrainAvailable)
                    {
                        string err = "You are using absolute coordinates but flat terrain  - ESC = Exit";
                        Console.WriteLine(err);
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
                }
                else
                {
                    Program.PS_Absolute_Height[i] = false;
                }

                Program.PS_Z[i] = (float)Math.Abs(PQ[i].Z1);
                
                if (i <= 5)//20260202
                {
                    Console.WriteLine($"[DBG] PS_Z[{i}]={Program.PS_Z[i]}");
                }

                Program.PS_V[i] = PQ[i].V;
                Program.PS_VDirX[i] = PQ[i].VDirX;//20260408 point source velocity direction x
                Program.PS_VDirY[i] = PQ[i].VDirY;//20260408 point source velocity direction y
                Program.PS_VDirZ[i] = PQ[i].VDirZ;//20260408 point source velocity direction z
                Program.PS_D[i] = PQ[i].D;
                Program.PS_T[i] = PQ[i].T;
                Program.PS_SG[i] = (byte)PQ[i].SG;
                Program.PS_V_Dep[i] = PQ[i].Vdep;
                Program.PS_V_sed[i] = PQ[i].Vsed;
                Program.PS_Mode[i] = PQ[i].Mode;
                Program.PS_ER_Dep[i] = (float)(PQ[i].ER_dep);
                Program.PS_TimeSeriesTemperature[i] = PQ[i].TimeSeriesTemperature;
                Program.PS_TimeSeriesVelocity[i] = PQ[i].TimeSeriesVelocity;
            }

            string info = "Total number of point sources: " + countrealsources.ToString();
            Console.WriteLine(info);
            ProgramWriters.LogfileGralCoreWrite(info);

            string unit = "[kg/h]: ";
            if (Program.Odour == true)
            {
                unit = "[MOU/h]: ";
            }

            info = "Total emission " + unit + (totalemission).ToString("0.000");
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

            Program.PS_effqu = new float[Program.PS_Count + 1];

            PQ.Clear();
            PQ.TrimExcess();
            Dep = null;
        }
    }
}
