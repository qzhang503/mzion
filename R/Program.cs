using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;


namespace MYSPACE
{
    public static class StringArrayExtensions
    {
        public static int NthIndexOf(this string[] array, string value, int n)
        {
            int count = 0;
            for (int i = 0; i < array.Length; i++)
            {
                if (array[i] == value)
                {
                    count++;
                    if (count == n)
                    {
                        return i;
                    }
                }
            }
            return -1;
        }
    }

    public static class IRawDataPlusExtension
    {
        public static void WriteSpectrum(this IRawDataPlus rawFile, string outname, List<int> L)
        {
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(outname))
            {
                double isoMass = 0.0, isoWidth = 0.0;
                string msOrder = "0";

                file.WriteLine("RAW\n{0}\n", Path.GetFileName(rawFile.FileName));

                foreach (int scanNumber in L)
                {
                    var scanStatistics = rawFile.GetScanStatsForScanNumber(scanNumber);
                    var scanEvent = rawFile.GetScanEventForScanNumber(scanNumber);

                    file.WriteLine("SCAN\n{0}", scanNumber);
                    file.WriteLine("TIME\n{0}", Math.Round(Convert.ToDouble(scanStatistics.StartTime), 6));

                    msOrder = scanEvent.MSOrder.ToString().ToLower();

                    if (msOrder == "ms2")
                    {
                        msOrder = "2";
                    }
                    else if (msOrder == "ms")
                    {
                        msOrder = "1";
                    }
                    else
                    {
                        msOrder = "0";
                    }

                    file.WriteLine("MSORDER\n{0}", msOrder);

                    if (msOrder == "1")
                    {
                        isoMass = 0;
                        isoWidth = 0;
                    }
                    else
                    {
                        var reaction0 = scanEvent.GetReaction(0); // not for PROFILE MS1
                        isoMass = Math.Round(Convert.ToDouble(reaction0.PrecursorMass), 4);
                        isoWidth = Math.Round(Convert.ToDouble(reaction0.IsolationWidth), 4);
                    }

                    var centroidStream = rawFile.GetCentroidStream(scanNumber, false);
                    var NcentroidPeaks = centroidStream.Length;

                    if (NcentroidPeaks > 0)
                    {
                        double[] masses = new double[NcentroidPeaks];
                        double[] intensities = new double[NcentroidPeaks];

                        for (int i = 0; i < NcentroidPeaks; i++)
                        {
                            masses[i] = Math.Round(Convert.ToDouble(centroidStream.Masses[i]), 4);
                            intensities[i] = Math.Round(Convert.ToDouble(centroidStream.Intensities[i]), 1);
                        }

                        file.WriteLine("ISOMASS\n{0}", isoMass);
                        file.WriteLine("ISOWIDTH\n{0}", isoWidth);
                        file.WriteLine("NPEAKS\n{0}", NcentroidPeaks);
                        file.WriteLine("X\n{0}", string.Join(" ", masses));
                        file.WriteLine("Y\n{0}", string.Join(" ", intensities));
                        file.WriteLine("TITLE\nFile: {0}; scans: {1}\n", rawFile.FileName, scanNumber);
                    }
                    else
                    {
                        file.WriteLine("ISOMASS\n0");
                        file.WriteLine("ISOWIDTH\n0");
                        file.WriteLine("NPEAKS\n0");
                        file.WriteLine("X\n0");
                        file.WriteLine("Y\n0");
                        file.WriteLine("TITLE\nFile: {0}; scans: {1}\n", rawFile.FileName, scanNumber);
                    }

                }
            }

            return;
        }

    }
}


namespace MyApp
{
    using MYSPACE;

    internal static class Program
    {
        private static void Main(string[] args)
        {
            /*
             #if DEBUG
                        args = new[] { "mypath/filename.raw",
                        "~/AppData/Local/Temp/foo/foo.peaks"
                        };
            #endif 
            */

            string filename = args[0];
            string outname  = args[1];

            if (string.IsNullOrEmpty(filename))
            {
                Console.WriteLine("No RAW file specified!");
                return;
            }

            if (string.IsNullOrEmpty(outname))
            {
                Console.WriteLine("No output file specified!");
                return;
            }

            try
            {
                var rawFile = RawFileReaderAdapter.FileFactory(filename);

                if (!rawFile.IsOpen || rawFile.IsError)
                {
                    Console.WriteLine("Unable to access the RAW file using the RawFileReader class!");
                    return;
                }

                if (rawFile.IsError)
                {
                    return;
                }

                if (rawFile.InAcquisition)
                {
                    return;
                }

                rawFile.SelectInstrument(Device.MS, 1);

                int firstScanNumber = rawFile.RunHeaderEx.FirstSpectrum;
                int lastScanNumber = rawFile.RunHeaderEx.LastSpectrum;
                List<int> scans = Enumerable.Range(firstScanNumber, lastScanNumber).ToList();

                if (scans.Count > 0)
                {
                    rawFile.WriteSpectrum(outname, scans);
                }

                return;
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error accessing RAWFileReader library! - " + ex.Message);
            }
        }

    }
}

