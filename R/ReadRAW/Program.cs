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
        public static void WriteSpectrum(this IRawDataPlus rawFile, string filename, List<int> L)
        {
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(filename))
            {
                int IndexMonoMZ = 8, IndexIsoWidth = 12;
                double isoMass = 0.0, isoWidth = 0.0;
                string msOrder = "0";
                bool isCentroidMS2 = true;

                file.WriteLine("RAW\n{0}\n", Path.GetFileName(rawFile.FileName));

                foreach (int scanNumber in L)
                {
                    var scanStatistics = rawFile.GetScanStatsForScanNumber(scanNumber);
                    var scanEvent = rawFile.GetScanEventForScanNumber(scanNumber);

                    file.WriteLine("SCAN\n{0}", scanNumber);
                    file.WriteLine("TIME\n{0}", Math.Round(Convert.ToDouble(scanStatistics.StartTime), 6));

                    // MS order
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

                    // Centroid MS2
                    isCentroidMS2 = scanStatistics.IsCentroidScan && msOrder == "2";

                    if (isCentroidMS2)
                    {
                        var reaction0 = scanEvent.GetReaction(0); // only available with centroid scans
                        isoMass = Math.Round(Convert.ToDouble(reaction0.PrecursorMass), 4);
                        isoWidth = Math.Round(Convert.ToDouble(reaction0.IsolationWidth), 4);
                    }

                    var centroidStream = rawFile.GetCentroidStream(scanNumber, false);
                    var NcentroidPeaks = centroidStream.Length;

                    if (NcentroidPeaks > 0)
                    {
                        var scanTrailer = rawFile.GetTrailerExtraInformation(scanNumber);
                        IndexMonoMZ = scanTrailer.Labels.NthIndexOf("Monoisotopic M/Z:", 1);
                        IndexIsoWidth = scanTrailer.Labels.NthIndexOf("MS2 Isolation Width:", 1);

                        double[] masses = new double[NcentroidPeaks];
                        double[] intensities = new double[NcentroidPeaks];

                        for (int i = 0; i < NcentroidPeaks; i++)
                        {
                            masses[i] = Math.Round(Convert.ToDouble(centroidStream.Masses[i]), 4);
                            intensities[i] = Math.Round(Convert.ToDouble(centroidStream.Intensities[i]), 1);
                        }

                        if (isCentroidMS2)
                        {
                            file.WriteLine("ISOMASS\n{0}", isoMass);
                            file.WriteLine("ISOWIDTH\n{0}", isoWidth);
                        }
                        else
                        {
                            if (IndexMonoMZ >= 0)
                            {
                                file.WriteLine("ISOMASS\n{0}", scanTrailer.Values[IndexMonoMZ]); // "-1.00" for DIA
                            }
                            else
                            {
                                file.WriteLine("ISOMASS\n0");
                            }

                            if (IndexIsoWidth >= 0)
                            {
                                file.WriteLine("ISOWIDTH\n{0}", scanTrailer.Values[IndexIsoWidth]); // "-1.00" for DIA
                            }
                            else
                            {
                                file.WriteLine("ISOWIDTH\n0");
                            }
                        }
                        
                        file.WriteLine("NPEAKS\n{0}", NcentroidPeaks);
                        file.WriteLine("X\n{0}", string.Join(",", masses));
                        file.WriteLine("Y\n{0}", string.Join(",", intensities));
                        // file.WriteLine("LOWMASS\n{0}", masses[0]);
                        // file.WriteLine("HIGHMASS\n{0}", masses[masses.Length - 1]);
                        file.WriteLine("TITLE\nFile: {0}; scans: {1}\n", rawFile.FileName, scanNumber);
                    }
                    else
                    {
                        file.WriteLine("ISOMASS\n0");
                        file.WriteLine("ISOWIDTH\n0");
                        file.WriteLine("NPEAKS\n0");
                        file.WriteLine("X\n0");
                        file.WriteLine("Y\n0");
                        // file.WriteLine("LOWMASS\n0");
                        // file.WriteLine("HIGHMASS\n0");
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

            string filename = string.Empty;

            filename = args[0];

            if (string.IsNullOrEmpty(filename))
            {
                Console.WriteLine("No RAW file specified!");
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
                    rawFile.WriteSpectrum(args[1], scans);
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

