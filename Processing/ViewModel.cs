using Microsoft.Win32;
using System;
using System.IO;
using System.Text;
using System.Windows;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace Processing
{
    sealed class ViewModel : NotificationObject
    {
        #region edge detect property
        private double xmin, xmax, ymin, ymax;
        private double[,] anomalyData, resultData;

        private bool isFrqDomain;

        public bool IsFrqDomain
        {
            get { return isFrqDomain; }
            set
            {
                isFrqDomain = value;
                RaisePropertyChanged("IsFrqDomain");
            }
        }

        private int winSize;

        public int WinSize
        {
            get { return winSize; }
            set
            {
                winSize = value;
                RaisePropertyChanged("WinSize");
            }
        }

        private double theta;

        public double Theta
        {
            get { return theta; }
            set
            {
                theta = value;
                RaisePropertyChanged("Theta");
            }
        }

        private double alpha;

        public double Alpha
        {
            get { return alpha; }
            set
            {
                alpha = value;
                RaisePropertyChanged("Alpha");
            }
        }

        private Complex2D frqData;

        public Complex2D FrqData
        {
            get
            {
                if (frqData == null)
                    frqData = FFT.fft_2D(new Complex2D(EdgeDetection.MinCurveExpand(anomalyData)));
                return frqData;
            }
            set { frqData = value; }
        }

        public double[,] Tempxdrdata
        {
            get
            {
                if (isFrqDomain)
                    return XdrfData;
                else
                    return XdrData;
            }
        }

        public double[,] Tempydrdata
        {
            get
            {
                if (isFrqDomain)
                    return YdrfData;
                else
                    return YdrData;
            }
        }

        public double[,] Tempzdrdata
        {
            get
            {
                if (isFrqDomain)
                    return ZdrfData;
                else
                    return ZdrData;
            }
        }

        private double[,] xdrData;

        public double[,] XdrData
        {
            get
            {
                if (xdrData == null)
                    xdrData = EdgeDetection.Xdr(anomalyData);
                return xdrData;
            }
            set { xdrData = value; }
        }

        private double[,] ydrData;

        public double[,] YdrData
        {
            get
            {
                if (ydrData == null)
                    ydrData = EdgeDetection.Ydr(anomalyData);
                return ydrData;
            }
            set { ydrData = value; }
        }

        private double[,] zdrData;

        public double[,] ZdrData
        {
            get
            {
                if (zdrData == null)
                    zdrData = EdgeDetection.Zdr(anomalyData);
                return zdrData;
            }
            set { zdrData = value; }
        }

        private double[,] xdrfData;

        public double[,] XdrfData
        {
            get
            {
                if (xdrfData == null)
                    xdrfData = EdgeDetection.ExpandRestore(FFT.idft_2D(EdgeDetection.Xdrf(FrqData)).GetRe());
                return xdrfData;
            }
            set { xdrfData = value; }
        }

        private double[,] ydrfData;

        public double[,] YdrfData
        {
            get
            {
                if (ydrfData == null)
                    ydrfData = EdgeDetection.ExpandRestore(FFT.idft_2D(EdgeDetection.Ydrf(FrqData)).GetRe());
                return ydrfData;
            }
            set { ydrfData = value; }
        }

        private double[,] zdrfData;

        public double[,] ZdrfData
        {
            get
            {
                if (zdrfData == null)
                    zdrfData = EdgeDetection.ExpandRestore(FFT.idft_2D(EdgeDetection.Zdrf(FrqData)).GetRe());
                return zdrfData;
            }
            set { zdrfData = value; }
        }

        private double[,] vxxData;

        public double[,] VxxData
        {
            get
            {
                if (vxxData == null)
                    vxxData = EdgeDetection.ExpandRestore(FFT.idft_2D(EdgeDetection.Vxx(FrqData)).GetRe());
                return vxxData;
            }
            set { vxxData = value; }
        }

        private double[,] vyyData;

        public double[,] VyyData
        {
            get
            {
                if (vyyData == null)
                    vyyData = EdgeDetection.ExpandRestore(FFT.idft_2D(EdgeDetection.Vyy(FrqData)).GetRe());
                return vyyData;
            }
            set { vyyData = value; }
        }

        private double[,] vxyData;

        public double[,] VxyData
        {
            get
            {
                if (vxyData == null)
                    vxyData = EdgeDetection.ExpandRestore(FFT.idft_2D(EdgeDetection.Vxy(FrqData)).GetRe());
                return vxyData;
            }
            set { vxyData = value; }
        }

        public DelegateCommand OpenFileCommand { get; set; }
        public DelegateCommand SaveFileCommand { get; set; }
        public DelegateCommand EdgeDetectCommand { get; set; }
        #endregion

        #region anomaly separation property.

        private double[,] secondData, secondResultData;

        private int power;

        public int Power
        {
            get { return power; }
            set
            {
                power = value;
                RaisePropertyChanged("Power");
            }
        }

        private int num;

        public int Num
        {
            get { return num; }
            set
            {
                num = value;
                RaisePropertyChanged("Num");
            }
        }


        private int radius;

        public int Radius
        {
            get { return radius; }
            set
            {
                radius = value;
                RaisePropertyChanged("Radius");
            }
        }

        private int radiusMin;

        public int RadiusMin
        {
            get { return radiusMin; }
            set
            {
                radiusMin = value;
                RaisePropertyChanged("RadiusMin");
            }
        }

        private int radiusMax;

        public int RadiusMax
        {
            get { return radiusMax; }
            set
            {
                radiusMax = value;
                RaisePropertyChanged("RadiusMax");
            }
        }

        private int radiusTick;

        public int RadiusTick
        {
            get { return radiusTick; }
            set
            {
                radiusTick = value;
                RaisePropertyChanged("RadiusTick");
            }
        }

        private int times;

        public int Times
        {
            get { return times; }
            set
            {
                times = value;
                RaisePropertyChanged("Times");
            }
        }

        private double project;

        public double Project
        {
            get { return project; }
            set
            {
                project = value;
                RaisePropertyChanged("Project");
            }
        }

        private int length;

        public int Length
        {
            get { return length; }
            set
            {
                length = value;
                RaisePropertyChanged("Length");
            }
        }

        private int expendMethod;

        public int ExpendMethod
        {
            get { return expendMethod; }
            set
            {
                expendMethod = value;
                RaisePropertyChanged("ExpendMethod");
            }
        }

        private double yantuogaodu;

        public double 延拓高度
        {
            get { return yantuogaodu; }
            set
            {
                yantuogaodu = value;
                RaisePropertyChanged("延拓高度");
            }
        }


        public double[,] FiltData
        {
            get
            {
                if (expendMethod == 0)
                {
                    return Separation.Filter(anomalyData, 0, project, length);
                }
                else if (expendMethod == 1)
                {
                    return Separation.Filter(anomalyData, 1, project, length);
                }
                else
                {
                    return Separation.Filter(anomalyData, ZdrfData, project, length);
                }
            }
        }


        public DelegateCommand SeparateCommand { get; set; }
        #endregion

        #region render
        private Colormap colormap;

        private double anomalyMin;

        public double AnomalyMin
        {
            get { return anomalyMin; }
            set
            {
                anomalyMin = value;
                RaisePropertyChanged("AnomalyMin");
            }
        }

        private double anomalyMax;

        public double AnomalyMax
        {
            get { return anomalyMax; }
            set
            {
                anomalyMax = value;
                RaisePropertyChanged("AnomalyMax");
            }
        }

        private WriteableBitmap anomalyMap;

        public WriteableBitmap AnomalyMap
        {
            get { return anomalyMap; }
            set
            {
                anomalyMap = value;
                RaisePropertyChanged("AnomalyMap");
            }
        }

        private double secondMin;

        public double SecondMin
        {
            get { return secondMin; }
            set
            {
                secondMin = value;
                RaisePropertyChanged("SecondMin");
            }
        }

        private double secondMax;

        public double SecondMax
        {
            get { return secondMax; }
            set
            {
                secondMax = value;
                RaisePropertyChanged("SecondMax");
            }
        }

        private WriteableBitmap secondMap;

        public WriteableBitmap SedondMap
        {
            get { return secondMap; }
            set
            {
                secondMap = value;
                RaisePropertyChanged("SedondMap");
            }
        }

        private double resultMin;

        public double ResultMin
        {
            get { return resultMin; }
            set
            {
                resultMin = value;
                RaisePropertyChanged("ResultMin");
            }
        }

        private double resutlMax;

        public double ResultMax
        {
            get { return resutlMax; }
            set
            {
                resutlMax = value;
                RaisePropertyChanged("ResultMax");
            }
        }

        private WriteableBitmap resultMap;

        public WriteableBitmap ResultMap
        {
            get { return resultMap; }
            set
            {
                resultMap = value;
                RaisePropertyChanged("ResultMap");
            }
        }

        private double secondResultMin;

        public double SecondResultMin
        {
            get { return secondResultMin; }
            set
            {
                secondResultMin = value;
                RaisePropertyChanged("SecondResultMin");
            }
        }

        private double secondResultMax;

        public double SecondResultMax
        {
            get { return secondResultMax; }
            set
            {
                secondResultMax = value;
                RaisePropertyChanged("SecondResultMax");
            }
        }

        private WriteableBitmap secondResultMap;

        public WriteableBitmap SecondResultMap
        {
            get { return secondResultMap; }
            set
            {
                secondResultMap = value;
                RaisePropertyChanged("SecondResultMap");
            }
        }
        #endregion

        public DelegateCommand GridMinusCommand { get; set; }
        public DelegateCommand FillRestoreCommand { get; set; }

        public ViewModel()
        {
            WinSize = 5;
            Power = 3;
            RadiusMin = 2;
            radiusMax = 20;
            RadiusTick = 2;
            Times = 3;
            project = 0.5;
            length = 50;
            OpenFileCommand = new DelegateCommand(new Action<object>(OpenFile));
            SaveFileCommand = new DelegateCommand(new Action<object>(SaveFile));
            EdgeDetectCommand = new DelegateCommand(new Action<object>(EdgeDetect));
            SeparateCommand = new DelegateCommand(new Action<object>(Separate));
            GridMinusCommand = new DelegateCommand(new Action<object>(GridMinus));
            FillRestoreCommand = new DelegateCommand(new Action<object>(FillRestore));

            colormap = new Colormap();
            colormap.SetColors();
        }

        private void OpenFile(object parameter)
        {
            var ofd = new OpenFileDialog
            {
                Filter = "Grid|*.grd"
            };

            if (ofd.ShowDialog() == true)
            {
                int rows, cols;
                double xtick, ytick, zmin, zmax;
                var para = Convert.ToString(parameter);
                if (para == "Anomaly")
                {
                    anomalyData = OpenFile(ofd.OpenFile(), out rows, out cols, out xtick, out ytick, out zmin, out zmax);

                    AnomalyMin = zmin;
                    AnomalyMax = zmax;
                    AnomalyMap = GetMap(anomalyData, zmin, zmax);

                    EdgeDetection.Initialize(cols, rows, xtick, ytick, anomalyData);

                    frqData = null;
                    xdrData = null;
                    ydrData = null;
                    zdrData = null;
                    xdrfData = null;
                    ydrfData = null;
                    zdrfData = null;
                    vxxData = null;
                    vyyData = null;
                    vxyData = null;
                }
                else if (para == "Local")
                {
                    secondData = OpenFile(ofd.OpenFile(), out rows, out cols, out xtick, out ytick, out zmin, out zmax);

                    SecondMin = zmin;
                    SecondMax = zmax;
                    SedondMap = GetMap(secondData, zmin, zmax);
                }
            }

        }

        private WriteableBitmap GetMap(double[,] data, double zmin, double zmax)
        {
            int rows = data.GetLength(0);
            int cols = data.GetLength(1);
            var map = new WriteableBitmap(cols, rows, 96, 96, PixelFormats.Rgb24, null);
            var pixels = new byte[map.PixelHeight, map.BackBufferStride];
            int index = 0;
            for (int i = 0; i < rows; i++)
            {
                index = 0;
                for (int j = 0; j < cols; j++)
                {
                    var color = colormap.GetColor(data[i, j], zmin, zmax);
                    pixels[i, index++] = color.R;
                    pixels[i, index++] = color.G;
                    pixels[i, index++] = color.B;
                }
            }
            map.WritePixels(new Int32Rect(0, 0, cols, rows), pixels, map.BackBufferStride, 0);
            return map;
        }

        private double[,] OpenFile(Stream stream, out int rows, out int cols, out double xtick, out double ytick, out double zmin, out double zmax)
        {
            double[,] data;

            using (var sr = new StreamReader(stream))
            {
                string[] strs;
                char[] ch = new char[] { ' ' };

                sr.ReadLine();//msaa
                strs = sr.ReadLine().Split(ch);
                rows = int.Parse(strs[0]);
                cols = int.Parse(strs[1]);
                data = new double[rows, cols];

                strs = sr.ReadLine().Split(ch);
                xmin = double.Parse(strs[0]);
                xmax = double.Parse(strs[1]);
                xtick = (xmax - xmin) / (cols - 1);

                strs = sr.ReadLine().Split(ch);
                ymin = double.Parse(strs[0]);
                ymax = double.Parse(strs[1]);
                ytick = (ymax - ymin) / (rows - 1);

                strs = sr.ReadLine().Split(ch);
                zmin = double.Parse(strs[0]);
                zmax = double.Parse(strs[1]);

                for (int i = 0; i < rows; i++)
                {
                    int j = 0;
                    while (j < cols)
                    {
                        strs = sr.ReadLine().Split(ch, StringSplitOptions.RemoveEmptyEntries);
                        for (int col = 0; col < strs.Length; col++)
                        {
                            data[i, j] = double.Parse(strs[col]);
                            j++;
                        }
                    }
                }
            }

            return data;
        }

        private void SaveFile(object parameter)
        {
            var sfd = new SaveFileDialog()
            {
                Filter = "Grid|*.grd",
            };

            if (sfd.ShowDialog() == true)
            {
                var para = Convert.ToString(parameter);
                using (var fileStream = File.Create(sfd.FileName))
                {
                    if (para == "Result")
                    {
                        SaveFile(fileStream, resultData);
                    }
                    else if (para == "Local")
                    {
                        SaveFile(fileStream, secondResultData);
                    }
                }

            }
        }

        private void SaveFile(FileStream fileStream, double[,] data)
        {
            int rows = data.GetLength(0);
            int cols = data.GetLength(1);
            using (var sw = new StreamWriter(fileStream))
            {
                sw.WriteLine("DSAA");
                sw.WriteLine(rows.ToString() + " " + cols.ToString());
                sw.WriteLine(xmin.ToString() + " " + xmax.ToString());
                sw.WriteLine(ymin.ToString() + " " + ymax.ToString());
                var zmin = EdgeDetection.GetDataMin(data);
                var zmax = EdgeDetection.GetDataMax(data);
                sw.WriteLine(zmin.ToString() + " " + zmax.ToString());

                var line = new StringBuilder();
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        line.Append(data[i, j]);
                        line.Append(" ");
                    }

                    sw.WriteLine(line.ToString());
                    line.Clear();
                }
            }
        }

        private void GridMinus(object parameter)
        {
            int rows = anomalyData.GetLength(0);
            int cols = anomalyData.GetLength(1);
            resultData = new double[rows, cols];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    resultData[i, j] = anomalyData[i, j] - secondData[i, j];
                }
            }
            var zmin = EdgeDetection.GetDataMin(resultData);
            var zmax = EdgeDetection.GetDataMax(resultData);
            ResultMin = zmin;
            ResultMax = zmax;
            ResultMap = GetMap(resultData, zmin, zmax);
        }

        private void GridSimplify(object parameter)
        {
            int rows = anomalyData.GetLength(0);
            int cols = anomalyData.GetLength(1);
            resultData = new double[51, 51];
            for (int i = 0; i < rows; i += 4)
            {
                for (int j = 0; j < cols; j += 4)
                {
                    resultData[i / 4, j / 4] = anomalyData[i, j];
                }
            }
            var zmin = EdgeDetection.GetDataMin(resultData);
            var zmax = EdgeDetection.GetDataMax(resultData);
            ResultMin = zmin;
            ResultMax = zmax;
            ResultMap = GetMap(resultData, zmin, zmax);
        }

        private void EdgeDetect(object parameter)
        {
            int index = Convert.ToInt32(parameter);
            if (index == 0)
                resultData = Tempxdrdata;
            else if (index == 1)
                resultData = Tempydrdata;
            else if (index == 2)
                resultData = Tempzdrdata;
            else if (index == 3)
                resultData = EdgeDetection.Thdr(Tempxdrdata, Tempydrdata);
            else if (index == 4)
                resultData = EdgeDetection.Asm(Tempxdrdata, Tempydrdata, Tempzdrdata);
            else if (index == 5)
                resultData = EdgeDetection.Tilt(Tempzdrdata, EdgeDetection.Thdr(Tempxdrdata, Tempydrdata));
            else if (index == 6)
                resultData = EdgeDetection.Theta(EdgeDetection.Thdr(Tempxdrdata, Tempydrdata), EdgeDetection.Asm(Tempxdrdata, Tempydrdata, Tempzdrdata));
            else if (index == 7)
                resultData = EdgeDetection.Nstd(Tempxdrdata, Tempydrdata, Tempzdrdata, winSize);
            else if (index == 8)
                resultData = EdgeDetection.Gdo(Tempxdrdata, Tempydrdata, Tempzdrdata, alpha, theta);
            else if (index == 9)
                resultData = EdgeDetection.Ita2(Tempxdrdata, Tempydrdata, Tempzdrdata);
            else if (index == 10)
                resultData = EdgeDetection.ETA1(Tempxdrdata, Tempydrdata, Tempzdrdata);
            else if (index == 11)
                resultData = EdgeDetection.HFTilt(Tempxdrdata, Tempydrdata, VxxData, VyyData, VxyData);
            else if (index == 12)
                resultData = EdgeDetection.HFTheta(Tempxdrdata, Tempydrdata, VxxData, VyyData, VxyData);
            else if (index == 13)
            {
                var asmData = EdgeDetection.Asm(Tempxdrdata, Tempydrdata, Tempzdrdata);
                double[,] asmZdrData, asmXdrData, asmYdrData;
                if (isFrqDomain)
                {
                    var asmFrqData = FFT.fft_2D(new Complex2D(EdgeDetection.MinCurveExpand(asmData)));
                    asmZdrData = EdgeDetection.ExpandRestore(FFT.idft_2D(EdgeDetection.Zdrf(asmFrqData)).GetRe());
                    asmYdrData = EdgeDetection.ExpandRestore(FFT.idft_2D(EdgeDetection.Ydrf(asmFrqData)).GetRe());
                    asmXdrData = EdgeDetection.ExpandRestore(FFT.idft_2D(EdgeDetection.Xdrf(asmFrqData)).GetRe());
                }
                else
                {
                    asmZdrData = EdgeDetection.Zdr(asmData);
                    asmYdrData = EdgeDetection.Ydr(asmData);
                    asmXdrData = EdgeDetection.Xdr(asmData);
                }
                var asmThdrData = EdgeDetection.Thdr(asmXdrData, asmYdrData);
                resultData = EdgeDetection.Tilt(asmZdrData, asmThdrData);
            }
            else if (index == 14)
            {
                var tiltData = EdgeDetection.Tilt(Tempzdrdata, EdgeDetection.Thdr(Tempxdrdata, Tempydrdata));
                double[,] tiltXdrData, tiltYdrData;
                if (isFrqDomain)
                {
                    var tiltFrqData = FFT.fft_2D(new Complex2D(EdgeDetection.MinCurveExpand(tiltData)));
                    tiltYdrData = EdgeDetection.ExpandRestore(FFT.idft_2D(EdgeDetection.Ydrf(tiltFrqData)).GetRe());
                    tiltXdrData = EdgeDetection.ExpandRestore(FFT.idft_2D(EdgeDetection.Xdrf(tiltFrqData)).GetRe());

                }
                else
                {
                    tiltYdrData = EdgeDetection.Ydr(tiltData);
                    tiltXdrData = EdgeDetection.Xdr(tiltData);
                }

                resultData = EdgeDetection.Thdr(tiltXdrData, tiltYdrData);
            }
            else if (index == 15)
                resultData = EdgeDetection.GuassNoise(anomalyData);
            else if (index == 16)
                resultData = EdgeDetection.Idct(EdgeDetection.Zdrf(EdgeDetection.Dct(anomalyData, false)));

            var zmin = EdgeDetection.GetDataMin(resultData);
            var zmax = EdgeDetection.GetDataMax(resultData);
            ResultMin = zmin;
            ResultMax = zmax;
            ResultMap = GetMap(resultData, zmin, zmax);
        }

        private void Separate(object parameter)
        {
            int rows = anomalyData.GetLength(0);
            int cols = anomalyData.GetLength(1);
            int index = Convert.ToInt32(parameter);

            if (index == 0)
                resultData = Separation.TrendAnalysis(anomalyData, power);
            else if (index == 1)
                resultData = Separation.FilterTrendAnalysis(anomalyData, FiltData, power);
            else if (index == 2)
                resultData = Separation.InterpolateCut(anomalyData, expendMethod, radius, times);
            else if (index == 3)
            {
                resultData = Separation.MoveAverage(anomalyData, secondData, expendMethod, radiusMin, radiusMax, radiusTick, out int bestRadius);
                Radius = bestRadius;
            }
            else if (index == 4)
                resultData = Separation.MoveAverage(anomalyData, expendMethod, radius);
            else if (index == 5)
                resultData = Separation.FilterMoveAverage(anomalyData, FiltData, expendMethod, radius);
            else if (index == 6)
                resultData = FiltData;
            else if (index == 7)
            {
                resultData = Separation.FilterProject(anomalyData, expendMethod, project);
            }
            else if (index == 8)
            {
                resultData = EdgeDetection.ExpandRestore(FFT.idft_2D(EdgeDetection.延拓(FrqData, 延拓高度)).GetRe());
            }

            var zmin = EdgeDetection.GetDataMin(resultData);
            var zmax = EdgeDetection.GetDataMax(resultData);
            ResultMin = zmin;
            ResultMax = zmax;
            ResultMap = GetMap(resultData, zmin, zmax);

            if (index != 6)
            {
                secondResultData = new double[rows, cols];
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        secondResultData[i, j] = anomalyData[i, j] - resultData[i, j];
                    }
                }
            }
            else
            {
                secondResultData = Separation.AverageFill(resultData);
            }

            zmin = EdgeDetection.GetDataMin(secondResultData);
            zmax = EdgeDetection.GetDataMax(secondResultData);
            SecondResultMin = zmin;
            SecondResultMax = zmax;
            SecondResultMap = GetMap(secondResultData, zmin, zmax);
        }

        private void FillRestore(object parameter)
        {
            int rows = anomalyData.GetLength(0);
            int cols = anomalyData.GetLength(1);

            resultData = new double[rows, cols];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    if (secondData[i, j] == 1.70141e38)
                    {
                        resultData[i, j] = 1.70141e38;
                    }
                    else
                    {
                        resultData[i, j] = anomalyData[i, j];
                    }
                }
            }

            var zmin = EdgeDetection.GetDataMin(resultData);
            var zmax = EdgeDetection.GetDataMax(resultData);
            ResultMin = zmin;
            ResultMax = zmax;
            ResultMap = GetMap(resultData, zmin, zmax);
        }
    }
}
