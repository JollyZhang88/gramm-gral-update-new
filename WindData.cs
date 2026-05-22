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


/// <summary>
/// Wind Data Class
/// </summary>
public class WindData
{
    public string Date { get; set; }
    public int IWet { get; set; }
    public int IDisp { get; set; }
    public int Stab_Class { get; set; }
    public string Hour { get; set; }
    public string Month { get; set; }
    public string Day { get; set; }
    // 20260101新增：秒级支持
    public int HourInt { get; set; }
    public int Minute { get; set; }
    public int Second { get; set; }
    public float HourFloat { get; set; }
    public float Vel { get; set; }
    public float Dir { get; set; }
    public int Frequency { get; set; }
}