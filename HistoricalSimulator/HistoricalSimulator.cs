/* Copyright (C) 2012 Fairmat SRL (info@fairmat.com, http://www.fairmat.com/)
 * Author(s): Francesco Biondi (francesco.biondi@fairmat.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using DVPLDOM;
using DVPLI;

namespace HistoricalSimulator
{
    /// <summary>
    /// Historical simulator allows to use historical trajectories of stocks as
    /// future realizations of underlying assets. In the current version trajectories
    /// are fetched from CSV or text files.
    /// </summary>
    /// <remarks>
    /// Futures versions will allow to retrieve trajectories from data providers and
    /// to simulate future realizations using bootstrapping techniques.
    /// </remarks>
    [Serializable]
    public class HistoricalSimulator : IFullSimulator, IExtensibleProcess, IParsable, IAssetSet
    {
        #region Fields
        /// <summary>
        /// The character used as separator for the CSV file.
        /// </summary>
        private static char[] elementSeparators = new char[] { ' ', ';', '\t' };

        /// <summary>
        /// A buffered copy of the lines of the file.
        /// </summary>
        [NonSerialized]
        private string[] fileLines;

        /// <summary>
        /// The starting index representing the entry point of the file.
        /// </summary>
        [NonSerialized]
        private int startIndex;

        /// <summary>
        /// Contains the mapping between the simulation dates and the line of the file to use for
        /// each date.
        /// </summary>
        [NonSerialized]
        private Dictionary<double, int> simulationDateIndexes;
        #endregion // Fields

        #region Properties
        /// <summary>
        /// Gets or sets the path to the CSV file.
        /// </summary>
        [PathSettingDescription("File path")]
        public string FilePath { get; set; }

        /// <summary>
        /// Gets or sets the start date used as entry point in the file.
        /// </summary>
        [SettingDescription("Start date")]
        public DateTime StartDate { get; set; }

        /// <summary>
        /// Gets or sets the operating mode: not used in the current version.
        /// </summary>
        [SettingDescription("Operating mode")]
        public OperatingMode OperatingMode { get; set; }

        #endregion // Properties

        #region Constructors
        /// <summary>
        /// Initializes the object.
        /// </summary>
        public HistoricalSimulator()
        {
            StartDate = DateTime.Now.Date;
        }
        #endregion // Constructors

        #region IFullSimulator implementation
        /// <summary>
        /// Simulate a realization of the stochastic process driven by the noise matrix noise
        /// this function is called once for realization.
        /// </summary>
        /// <param name="dates">
        /// The dates (in years fractions) at which the process must be simulated.
        /// </param>
        /// <param name="noise">The matrix of IID normal realizations.</param>
        /// <param name="outDynamic">Where the dynamic should be written.</param>
        public void Simulate(double[] dates, IReadOnlyMatrixSlice noise, IMatrixSlice outDynamic)
        {
            for (int i = 0; i < dates.Length; i++)
            {
                int dateIndex;
                if (!this.simulationDateIndexes.TryGetValue(dates[i], out dateIndex))
                    dateIndex = 0;

                string[] lineTokens = this.fileLines[dateIndex].Split(HistoricalSimulator.elementSeparators);
                for (int j = 0; j < lineTokens.Length - 1; j++)
                {
                    outDynamic[i, j] = DoubleHelper.FromString(lineTokens[j + 1]);
                }
            }
        }
        #endregion // IFullSimulator implementation

        #region IExtensibleProcess implementation
        /// <summary>
        /// Gets a value indicating whether FullSimulation is implemented, in this case it does
        /// so it always returns true.
        /// </summary>
        public bool ImplementsFullSimulation
        {
            get
            {
                return true;
            }
        }

        /// <summary>
        /// Gets a value indicating whether a Markov based simulation is implemented, in this
        /// case it doesn't so it always returns false.
        /// </summary>
        public bool ImplementsMarkovBasedSimulation
        {
            get
            {
                return false;
            }
        }

        /// <summary>
        /// Gets the ProcessInfo for this plug-in, in this case Historical Simulator.
        /// </summary>
        public ProcessInfo ProcessInfo
        {
            get
            {
                return new ProcessInfo("Historical Simulator");
            }
        }

        /// <summary>
        /// Called by Simulator after parse.
        /// Initializes here time-dependant but not state dependent variables.
        /// </summary>
        /// <param name="simulationDates">
        /// The dates at which the process realizations will be requested.
        /// </param>
        public void Setup(double[] simulationDates)
        {
            try
            {
                if (File.Exists(FilePath))
                {
                    // Store the lines of the file
                    this.fileLines = File.ReadAllLines(FilePath)
                                         .Where(s => s != string.Empty)
                                         .ToArray();

                    List<string> tmpList = new List<string>(this.fileLines);
                    tmpList.Sort((s1, s2) =>
                                 {
                                     string d1String = s1.Split(HistoricalSimulator.elementSeparators)[0];
                                     DateTime d1 = DateTime.Parse(d1String);
                                     string d2String = s2.Split(HistoricalSimulator.elementSeparators)[0];
                                     DateTime d2 = DateTime.Parse(d2String);

                                     return d1.CompareTo(d2);
                                 });
                    this.fileLines = tmpList.ToArray();

                    // Calculate the entry point in the file (the index of the starting line)
                    DateTime nearestDate = DateTime.Now.Date;
                    for (int i = 0; i < this.fileLines.Length; i++)
                    {
                        string[] lineTokens = this.fileLines[i].Split(HistoricalSimulator.elementSeparators);
                        DateTime date = DateTime.Parse(lineTokens[0]);
                        if (i == 0)
                            nearestDate = date;
                        else
                            nearestDate = NearestDate(StartDate, nearestDate, date);

                        if (nearestDate == date)
                            this.startIndex = i;
                    }

                    // Calculate the list of dates in the file
                    List<DateTime> dateList = new List<DateTime>();
                    for (int i = 0; i < this.fileLines.Length; i++)
                    {
                        string[] lineTokens = this.fileLines[i].Split(HistoricalSimulator.elementSeparators);
                        DateTime date = DateTime.Parse(lineTokens[0]);
                        dateList.Add(date);
                    }

                    // Calculate for each date the line of the file to use in order the get the
                    // values. Gets the nearest date between each simulation date and the dates
                    // present in the file (if the date is before the entry point use the entry
                    // point instead)
                    this.simulationDateIndexes = new Dictionary<double, int>();
                    foreach (double doubleDate in simulationDates)
                    {
                        DateTime date = Document.ActiveDocument.DefaultProject.GetDate(doubleDate);
                        TimeSpan dateDifference = date - Document.ActiveDocument.SimulationStartDate;
                        int dateIndex = NearestDateIndex(StartDate + dateDifference, dateList.ToArray());
                        this.simulationDateIndexes.Add(doubleDate, Math.Max(dateIndex, this.startIndex));
                    }
                }
                else
                    this.fileLines = null;
            }
            catch
            {
                this.fileLines = null;
            }
        }

        /// <summary>
        /// Gets the information required in order to allow the simulation to run.
        /// </summary>
        public SimulationInfo SimulationInfo
        {
            get
            {
                SimulationInfo simulationInfo = new SimulationInfo();
                simulationInfo.NoiseSize = 0;
                simulationInfo.LatentSize = 0;

                try
                {
                    // Use the number of columns of the first line of the file
                    // (At the moment don't check for an equal number of elements between each
                    // line)
                    if (File.Exists(FilePath))
                    {
                        string[] lines = File.ReadLines(FilePath).ToArray();
                        string[] lineTokens = lines[0].Split(HistoricalSimulator.elementSeparators);
                        simulationInfo.StateSize = lineTokens.Length - 1;
                    }
                }
                catch
                {
                    simulationInfo.StateSize = 0;
                }

                return simulationInfo;
            }
        }

        /// <summary>
        /// Creates a list of all the sub-objects that can be edited.
        /// </summary>
        /// <param name="recursive">The parameter is not used.</param>
        /// <returns>
        /// The created list with all the sub objects that can be edited (empty in this case).
        /// </returns>
        public List<IExportable> ExportObjects(bool recursive)
        {
            return new List<IExportable>();
        }
        #endregion // IExtensibleProcess implementation

        #region IParsable implementation
        /// <summary>
        /// Parses the process (in this case nothing has to be done).
        /// </summary>
        /// <param name="context">The project representing the context of the parsing.</param>
        /// <returns>true if the the parsing caused errors; otherwise false.</returns>
        public bool Parse(IProject context)
        {
            return false;
        }
        #endregion // IParsable implementation

        #region IAssetSet implementation
        /// <summary>
        /// Gets the list of Process indices belonging to the asset set.
        /// </summary>
        public int[] ProcessIndices
        {
            get
            {
                int[] indices = new int[SimulationInfo.StateSize];
                return indices;
            }
        }

        /// <summary>
        /// Gets the list of the indices in the state variables.
        /// </summary>
        public int[] StateIndices
        {
            get
            {
                int[] indices = new int[SimulationInfo.StateSize];
                return indices;
            }
        }
        #endregion // IAssetSet implementation

        /// <summary>
        /// Gets the nearest date to the reference date between the given dates.
        /// </summary>
        /// <param name="referenceDate">The date to use as reference.</param>
        /// <param name="dateA">The first date to check.</param>
        /// <param name="dateB">The second date to check.</param>
        /// <returns>The nearest date to the reference date.</returns>
        private DateTime NearestDate(DateTime referenceDate, DateTime dateA, DateTime dateB)
        {
            double absDifferenceA = Math.Abs((referenceDate - dateA).TotalDays);
            double absDifferenceB = Math.Abs((referenceDate - dateB).TotalDays);

            if (absDifferenceA < absDifferenceB)
                return dateA;
            else
                return dateB;
        }

        /// <summary>
        /// Gets the nearest date index to the reference date between the dates contained in the
        /// array.
        /// </summary>
        /// <param name="referenceDate">The date to use as reference.</param>
        /// <param name="dates">The array of dates to use.</param>
        /// <returns>The nearest date index to the reference date.</returns>
        private int NearestDateIndex(DateTime referenceDate, DateTime[] dates)
        {
            int retVal = 0;
            DateTime nearestDate = DateTime.Now.Date;
            for (int i = 0; i < dates.Length; i++)
            {
                if (i == 0)
                    nearestDate = dates[i];
                else
                    nearestDate = NearestDate(referenceDate, nearestDate, dates[i]);

                if (nearestDate == dates[i])
                    retVal = i;
            }

            return retVal;
        }
    }
}
