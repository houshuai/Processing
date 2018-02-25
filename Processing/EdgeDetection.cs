using System;
using System.Collections.Generic;

namespace Processing
{
    static class EdgeDetection
    {
        private static int ncols, nrows;
        private static double xtick, ytick;
        private static double[,] anomalyData;

        public static void Initialize(int _ncols, int _nrows, double _xtick, double _ytick, double[,] _anomalyData)
        {
            ncols = _ncols;
            nrows = _nrows;
            xtick = _xtick;
            ytick = _ytick;
            anomalyData = _anomalyData;
        }


        public static double[,] Xdr(double[,] data)
        {
            int nrows = data.GetLength(0);
            int ncols = data.GetLength(1);
            var ydrData = new double[nrows, ncols];
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 2; j < ncols - 2; j++)
                {
                    ydrData[i, j] = (data[i, j + 1] - data[i, j - 1]) / (10.0 * xtick / 1000.0) +
                        (data[i, j + 2] - data[i, j - 2]) / (5.0 * xtick / 1000.0);
                    //xderiData[i, j] = (data[i, j + 1] - data[i, j - 1]) / (2.0 * xtick / 1000.0);
                }
            }

            //最外层边界 向前、向后差分
            for (int i = 0; i < nrows; i++)
            {
                ydrData[i, 0] = (data[i, 1] - data[i, 0]) / (xtick / 1000.0);
                ydrData[i, ncols - 1] = (data[i, ncols - 1] - data[i, ncols - 2]) / (xtick / 1000.0);
            }

            //第二层边界 中心差分
            for (int i = 0; i < nrows; i++)
            {
                ydrData[i, 1] = (data[i, 2] - data[i, 0]) / (2.0 * xtick / 1000.0);
                ydrData[i, ncols - 2] = (data[i, ncols - 1] - data[i, ncols - 3]) / (2.0 * xtick / 1000.0);
            }
            return ydrData;
        }

        public static double[,] Ydr(double[,] data)
        {
            int nrows = data.GetLength(0);
            int ncols = data.GetLength(1);
            var xdrData = new double[nrows, ncols];
            for (int i = 0; i < ncols; i++)
            {
                for (int j = 2; j < nrows - 2; j++)
                {
                    xdrData[j, i] = (data[j + 1, i] - data[j - 1, i]) / (10.0 * ytick / 1000.0) +
                        (data[j + 2, i] - data[j - 2, i]) / (5.0 * ytick / 1000.0);
                    //xderiData[i, j] = (data[i, j + 1] - data[i, j - 1]) / (2.0 * xtick / 1000.0);
                }
            }

            //最外层边界 向前、向后差分
            for (int i = 0; i < ncols; i++)
            {
                xdrData[0, i] = (data[1, i] - data[0, i]) / (ytick / 1000.0);
                xdrData[nrows - 1, i] = (data[nrows - 1, i] - data[nrows - 2, i]) / (ytick / 1000.0);
            }

            //第二层边界 中心差分
            for (int i = 0; i < ncols; i++)
            {
                xdrData[1, i] = (data[2, i] - data[0, i]) / (2.0 * ytick / 1000.0);
                xdrData[nrows - 2, i] = (data[nrows - 1, i] - data[nrows - 3, i]) / (2.0 * ytick / 1000.0);
            }
            return xdrData;
        }

        public static double[,] Zdr(double[,] data)
        {
            int nrows = data.GetLength(0);
            int ncols = data.GetLength(1);
            var zdrData = new double[nrows, ncols];
            var R = new double[9] { 500, 1000, Math.Sqrt(2) * 1000, Math.Sqrt(5) * 1000, Math.Sqrt(8) * 1000,
                Math.Sqrt(13) * 1000, 5000, Math.Sqrt(50) * 1000, Math.Sqrt(136) * 1000 };

            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    //zdrData[i, j] = 10.0 * data[i, j] - (7.5 * getAve(i, j, 200) + 1.25 * getAve(i, j, 570) +
                    //    0.54 * getAve(i, j, 1060) + 0.26 * getAve(i, j, 1750) + 0.16 * getAve(i, j, 2730) +
                    //    0.09 * getAve(i, j, 4040) + 0.06 * getAve(i, j, 5630) + 0.04 * getAve(i, j, 7530) +
                    //    0.03 * getAve(i, j, 9990) + 0.02 * getAve(i, j, 13450) + 0.02 * getAve(i, j, 18360));
                    zdrData[i, j] = (1.937843 * data[i, j] - 0.464388 * GetAve(i, j, R[0]) -
                        0.518209 * GetAve(i, j, R[1]) - 0.274251 * GetAve(i, j, R[2]) -
                        0.293692 * GetAve(i, j, R[3]) - 0.06454 * GetAve(i, j, R[4]) -
                        0.079792 * GetAve(i, j, R[5]) - 0.058836 * GetAve(i, j, R[6]) -
                        0.078496 * GetAve(i, j, R[7]) - 0.020023 * GetAve(i, j, R[8]));
                }
            }
            return zdrData;
        }

        public static double GetAve(int y, int x, double r)
        {
            int xn = (int)(r / xtick);
            int yn = (int)(r / ytick);

            int xStart = x - xn >= 0 ? x - xn : 0;
            int xEnd = x + xn < ncols ? x + xn : ncols - 1;
            int yStart = y - yn >= 0 ? y - yn : 0;
            int yEnd = y + yn < nrows ? y + yn : nrows;

            double sum = 0;
            int count = 0;
            for (int i = xStart; i <= xEnd; i++)
            {
                double xlen = Math.Abs(i - x) * xtick;
                double ylen = Math.Sqrt(r * r - xlen * xlen);
                double rest = ylen % ytick;
                int n = (int)((ylen - rest) / ytick);
                if (y - n - 1 >= 0)
                {
                    sum += (anomalyData[y - n, i] * (ytick - rest) + anomalyData[y - n - 1, i] * rest) / ytick;
                    count += 1;
                }
                if (y + n + 1 < nrows)
                {
                    sum += (anomalyData[y + n, i] * (ytick - rest) + anomalyData[y + n + 1, i] * rest) / ytick;
                    count += 1;
                }

            }
            for (int i = yStart; i < yEnd; i++)
            {
                double ylen = Math.Abs(i - y) * ytick;
                double xlen = Math.Sqrt(r * r - ylen * ylen);
                double rest = xlen % xtick;
                int n = (int)((xlen - rest) / xtick);
                if (x - n - 1 >= 0)
                {
                    sum += (anomalyData[i, x - n] * (xtick - rest) + anomalyData[i, x - n - 1] * rest) / xtick;
                    count += 1;
                }
                if (x + n + 1 < ncols)
                {
                    sum += (anomalyData[i, x + n] * (xtick - rest) + anomalyData[i, x + n + 1] * rest) / xtick;
                    count += 1;
                }
            }

            return count == 0 ? 0 : sum / count;
        }

        public static Complex2D Xdrf(Complex2D frqData)
        {
            int height = frqData.Height;
            int whith = frqData.Width;
            var xdrfData = new Complex2D(height, whith);
            var coex = Coeff(whith, xtick);
            //var coey = Coeff(height, ytick);
            var im = new Complex(0, 1);
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < whith; j++)
                {
                    //var a = 2.0 * Math.PI * Math.Sqrt(coey[i] * coey[i] + coex[j] * coex[j]);
                    xdrfData[i, j] = coex[j] * frqData[i, j] * im;
                }
            }
            return xdrfData;
        }

        public static double[,] Xdrf(double[,] frqData)
        {
            int nrows = frqData.GetLength(0);
            int ncols = frqData.GetLength(1);
            var xdrfData = new double[nrows, ncols];
            var coex = Coeff(ncols, xtick);
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    ;
                    xdrfData[i, j] = coex[i] * frqData[i, j];
                }
            }
            return xdrfData;
        }

        public static Complex2D Ydrf(Complex2D frqData)
        {
            int height = frqData.Height;
            int whith = frqData.Width;
            var ydrfData = new Complex2D(height, whith);
            //var coex = Coeff(whith, xtick);
            var coey = Coeff(height, ytick);
            var im = new Complex(0, 1);
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < whith; j++)
                {
                    //var a = 2.0 * Math.PI * Math.Sqrt(coey[i] * coey[i] + coex[j] * coex[j]);
                    ydrfData[i, j] = coey[i] * frqData[i, j] * im;
                }
            }
            return ydrfData;
        }

        public static double[,] Ydrf(double[,] frqData)
        {
            int nrows = frqData.GetLength(0);
            int ncols = frqData.GetLength(1);
            var ydrfData = new double[nrows, ncols];
            var coey = Coeff(nrows, ytick);
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    ;
                    ydrfData[i, j] = coey[i] * frqData[i, j];
                }
            }
            return ydrfData;
        }

        public static Complex2D Zdrf(Complex2D frqData)
        {
            int height = frqData.Height;
            int whith = frqData.Width;
            var zdrfData = new Complex2D(height, whith);
            var coex = Coeff(whith, xtick);
            var coey = Coeff(height, ytick);
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < whith; j++)
                {
                    var a = Math.Sqrt(coey[i] * coey[i] + coex[j] * coex[j]);
                    zdrfData[i, j] = a * frqData[i, j];
                }
            }
            return zdrfData;
        }

        public static double[,] Zdrf(double[,] frqData)
        {
            int nrows = frqData.GetLength(0);
            int ncols = frqData.GetLength(1);
            var zdrfData = new double[nrows, ncols];
            var coex = Coeff(ncols, xtick);
            var coey = Coeff(nrows, ytick);
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    var a = Math.Sqrt(coey[i] * coey[i] + coex[j] * coex[j]);
                    zdrfData[i, j] = a * frqData[i, j];
                }
            }
            return zdrfData;
        }

        public static Complex2D Vxx(Complex2D frqData)
        {
            int height = frqData.Height;
            int whith = frqData.Width;
            var VxxData = new Complex2D(height, whith);
            var coex = Coeff(whith, xtick);
            var coey = Coeff(height, ytick);
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < whith; j++)
                {
                    VxxData[i, j] = -coey[i] * coey[i] * frqData[i, j] / Math.Sqrt(coey[i] * coey[i] + coex[j] * coex[j]);
                }
            }
            VxxData[0, 0] = new Complex(0, 0);
            VxxData[height / 2, 0] = new Complex(0, 0);
            VxxData[0, whith / 2] = new Complex(0, 0);
            VxxData[height / 2, whith / 2] = new Complex(0, 0);
            return VxxData;
        }

        public static Complex2D Vyy(Complex2D frqData)
        {
            int height = frqData.Height;
            int whith = frqData.Width;
            var VyyData = new Complex2D(height, whith);
            var coex = Coeff(whith, xtick);
            var coey = Coeff(height, ytick);
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < whith; j++)
                {
                    VyyData[i, j] = -coex[j] * coex[j] * frqData[i, j] / Math.Sqrt(coey[i] * coey[i] + coex[j] * coex[j]);
                }
            }
            VyyData[0, 0] = new Complex(0, 0);
            VyyData[height / 2, 0] = new Complex(0, 0);
            VyyData[0, whith / 2] = new Complex(0, 0);
            VyyData[height / 2, whith / 2] = new Complex(0, 0);
            return VyyData;
        }

        public static Complex2D Vxy(Complex2D frqData)
        {
            int height = frqData.Height;
            int whith = frqData.Width;
            var VxyData = new Complex2D(height, whith);
            var coex = Coeff(whith, xtick);
            var coey = Coeff(height, ytick);
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < whith; j++)
                {
                    VxyData[i, j] = -coey[i] * coex[j] * frqData[i, j] / Math.Sqrt(coey[i] * coey[i] + coex[j] * coex[j]);
                }
            }
            VxyData[0, 0] = new Complex(0, 0);
            VxyData[height / 2, 0] = new Complex(0, 0);
            VxyData[0, whith / 2] = new Complex(0, 0);
            VxyData[height / 2, whith / 2] = new Complex(0, 0);
            return VxyData;
        }

        private static double[] Coeff(int n, double tick)
        {
            var coe = new double[n];
            double f = 2.0 * Math.PI * 1000.0 / (n - 1) / tick;
            for (int i = 0; i <= n / 2; i++)
            {
                coe[i] = i * f;
                if (i != n / 2 && i != 0)
                    coe[n - i] = -i * f;
            }
            //for (int i = n / 2+1; i < n; i++)
            //{
            //    coe[i] = (i - n) * f;
            //}
            return coe;
        }


        public static double[,] Dct(double[,] data, bool is2PI)
        {
            int nrows = data.GetLength(0);
            int ncols = data.GetLength(1);
            var dctData = new double[nrows, ncols];

            if (is2PI) //水平导数
            {
                for (int m = 0; m < nrows; m++)
                {
                    for (int n = 0; n < ncols; n++)
                    {
                        dctData[m, n] = 0;
                        for (int i = 0; i < nrows; i++)
                        {
                            for (int j = 0; j < ncols; j++)
                            {
                                dctData[m, n] += data[i, j] * Math.Sin((2.0 * i + 1.0) * m * Math.PI / 2.0 / nrows) * Math.Sin((2.0 * j + 1.0) * n * Math.PI / 2.0 / ncols);
                            }
                        }
                        dctData[m, n] = dctData[m, n] * 2.0 / Math.Sqrt(nrows * ncols);
                    }
                }
            }
            else
            {
                for (int m = 0; m < nrows; m++)
                {
                    for (int n = 0; n < ncols; n++)
                    {
                        dctData[m, n] = 0;
                        for (int i = 0; i < nrows; i++)
                        {
                            for (int j = 0; j < ncols; j++)
                            {
                                dctData[m, n] += data[i, j] * Math.Cos((2.0 * i + 1.0) * m * Math.PI / 2.0 / nrows) * Math.Cos((2.0 * j + 1.0) * n * Math.PI / 2.0 / ncols);
                            }
                        }
                        dctData[m, n] = dctData[m, n] * 2.0 / Math.Sqrt(nrows * ncols);
                    }
                }
            }

            dctData[0, 0] /= 2.0;
            for (int n = 1; n < ncols; n++)
            {
                dctData[0, n] /= Math.Sqrt(2.0);
            }
            for (int m = 1; m < nrows; m++)
            {
                dctData[m, 0] /= Math.Sqrt(2.0);
            }

            return dctData;
        }

        public static double[,] Idct(double[,] frqData)
        {
            int nrows = frqData.GetLength(0);
            int ncols = frqData.GetLength(1);
            var data = new double[nrows, ncols];
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    data[i, j] = 0;
                    for (int m = 0; m < nrows; m++)
                    {
                        for (int n = 0; n < ncols; n++)
                        {
                            double c = 1.0;
                            if (m == 0)
                            {
                                if (n == 0)
                                    c = 1.0 / 2.0;
                                else
                                    c = 1.0 / Math.Sqrt(2.0);
                            }
                            else
                            {
                                if (n == 0)
                                    c = 1.0 / Math.Sqrt(2.0);
                            }
                            data[i, j] += c * frqData[m, n] * Math.Cos((2.0 * i + 1.0) * m * Math.PI / 2.0 / nrows) * Math.Cos((2.0 * j + 1.0) * n * Math.PI / 2.0 / ncols);
                        }
                    }
                    data[i, j] = data[i, j] * 2.0 / Math.Sqrt(nrows * ncols);
                }
            }
            return data;
        }

        public static double[,] Thdr(double[,] xdrData, double[,] ydrData)
        {
            int nrows = xdrData.GetLength(0);
            int ncols = xdrData.GetLength(1);
            var thdrData = new double[nrows, ncols];

            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    thdrData[i, j] = Math.Sqrt(xdrData[i, j] * xdrData[i, j] + ydrData[i, j] * ydrData[i, j]);
                }
            }

            return thdrData;
        }

        public static double[,] Asm(double[,] xdrData, double[,] ydrData, double[,] zdrData)
        {
            int nrows = xdrData.GetLength(0);
            int ncols = xdrData.GetLength(1);
            var asmData = new double[nrows, ncols];
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    asmData[i, j] = Math.Sqrt(Math.Pow(xdrData[i, j], 2) + Math.Pow(ydrData[i, j], 2) + Math.Pow(zdrData[i, j], 2));
                }
            }
            return asmData;
        }

        public static double[,] Tilt(double[,] zdrData, double[,] thdrData)
        {
            int nrows = zdrData.GetLength(0);
            int ncols = zdrData.GetLength(1);
            var tiltData = new double[nrows, ncols];
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    tiltData[i, j] = Math.Atan(zdrData[i, j] / thdrData[i, j]);
                }
            }
            return tiltData;
        }

        public static double[,] Theta(double[,] thdrData, double[,] asmData)
        {
            int nrows = thdrData.GetLength(0);
            int ncols = thdrData.GetLength(1);
            var thetaData = new double[nrows, ncols];
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    thetaData[i, j] = thdrData[i, j] / asmData[i, j];
                }
            }
            return thetaData;
        }

        public static double[,] Nstd(double[,] xdrData, double[,] ydrData, double[,] zdrData, int windowsize = 5)
        {
            int nrows = xdrData.GetLength(0);
            int ncols = xdrData.GetLength(1);
            var nstdData = new double[nrows, ncols];
            var xdrstdData = Std(xdrData, windowsize);
            var ydrstdData = Std(ydrData, windowsize);
            var zdrstdData = Std(zdrData, windowsize);
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    nstdData[i, j] = zdrstdData[i, j] / (xdrstdData[i, j] + ydrstdData[i, j] + zdrstdData[i, j]);
                }
            }

            return nstdData;
        }

        public static double[,] Std(double[,] data, int size)
        {
            int nrows = data.GetLength(0);
            int ncols = data.GetLength(1);
            var stdData = new double[nrows, ncols];
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    int xStart = j - size >= 0 ? j - size : 0;
                    int xEnd = j + size < ncols ? j + size : ncols - 1;
                    int yStart = i - size >= 0 ? i - size : 0;
                    int yEnd = i + size < nrows ? i + size : nrows - 1;
                    double miu = 0, sum = 0;
                    int n = 0;

                    for (int y = yStart; y < yEnd; y++)
                    {
                        for (int x = xStart; x < xEnd; x++)
                        {
                            miu += data[y, x];
                            n++;
                        }
                    }
                    miu /= n;
                    for (int y = yStart; y < yEnd; y++)
                    {
                        for (int x = xStart; x < xEnd; x++)
                        {
                            sum += Math.Pow((data[y, x] - miu), 2);
                        }
                    }
                    stdData[i, j] = Math.Sqrt(sum / n);
                }
            }

            return stdData;
        }

        public static double[,] Gdo(double[,] xdrData, double[,] ydrData, double[,] zdrData, double alpha = 0, double theta = 0)
        {
            int nrows = xdrData.GetLength(0);
            int ncols = xdrData.GetLength(1);
            var gdoData = new double[nrows, ncols];
            double sinAlpha = Math.Sin(alpha), sinTheta = Math.Sin(theta);
            double cosAlpha = Math.Cos(alpha), cosTheta = Math.Cos(theta);
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    gdoData[i, j] = ((xdrData[i, j] * sinAlpha + ydrData[i, j] * cosAlpha) * sinTheta + zdrData[i, j] * cosTheta) /
                        Math.Sqrt(Math.Pow(xdrData[i, j], 2) + Math.Pow(ydrData[i, j], 2) + Math.Pow(zdrData[i, j], 2));
                }
            }

            return gdoData;
        }

        public static double[,] Ita2(double[,] xdrData, double[,] ydrData, double[,] zdrData)
        {
            int nrows = xdrData.GetLength(0);
            int ncols = xdrData.GetLength(1);
            var ita2Data = new double[nrows, ncols];
            double[,] xdrxdrData, ydrydrData, zdrthdrData;
            xdrxdrData = Xdr(xdrData);
            ydrydrData = Ydr(ydrData);
            zdrthdrData = Thdr(Xdr(zdrData), Ydr(zdrData));

            var zdrzdrData = new double[nrows, ncols];
            double k;

            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    zdrzdrData[i, j] = -xdrxdrData[i, j] - ydrydrData[i, j];
                }
            }
            k = Math.Abs(GetDataMin(zdrData) / GetDataMax(zdrzdrData));

            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    ita2Data[i, j] = Math.Atan(k * zdrthdrData[i, j] / Math.Abs(zdrzdrData[i, j]));
                }
            }

            return ita2Data;
        }

        public static double[,] ETA1(double[,] xdrData, double[,] ydrData, double[,] zdrData)
        {
            int nrows = xdrData.GetLength(0);
            int ncols = xdrData.GetLength(1);
            var ETA1Data = new double[nrows, ncols];
            var zdrxdrData = Xdr(zdrData);
            var zdrydrData = Ydr(zdrData);
            var xdrxdrData = Xdr(xdrData);
            var ydrydrData = Ydr(ydrData);
            double c1 = Math.Abs(GetDataMin(zdrData) / GetDataMax(zdrData));
            var zdrzdrData = new double[nrows, ncols];

            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    zdrzdrData[i, j] = -xdrxdrData[i, j] - ydrydrData[i, j];
                }
            }

            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    ETA1Data[i, j] = Math.Atan(Math.Sqrt(zdrxdrData[i, j] * zdrxdrData[i, j] + zdrydrData[i, j] * zdrydrData[i, j] * c1) / Math.Abs(zdrzdrData[i, j]));
                }
            }

            return ETA1Data;
        }

        public static double[,] HFTilt(double[,] xdrData, double[,] ydrData, double[,] VxxData, double[,] VyyData, double[,] VxyData)
        {
            int nrows = xdrData.GetLength(0);
            int ncols = xdrData.GetLength(1);
            var HFTiltData = new double[nrows, ncols];
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    HFTiltData[i, j] = Math.Atan(Math.Sqrt(xdrData[i, j] * xdrData[i, j] + ydrData[i, j] * ydrData[i, j]) /
                        Math.Sqrt(VxxData[i, j] * VxxData[i, j] + VyyData[i, j] * VyyData[i, j] + 2.0 * VxyData[i, j] * VxyData[i, j]));
                }
            }
            return HFTiltData;
        }

        public static double[,] HFTheta(double[,] xdrData, double[,] ydrData, double[,] VxxData, double[,] VyyData, double[,] VxyData)
        {
            int nrows = xdrData.GetLength(0);
            int ncols = xdrData.GetLength(1);
            var HFThetaData = new double[nrows, ncols];
            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    HFThetaData[i, j] = Math.Sqrt(VxxData[i, j] * VxxData[i, j] + VyyData[i, j] * VyyData[i, j] + 2.0 * VxyData[i, j] * VxyData[i, j]) /
                        Math.Sqrt(xdrData[i, j] * xdrData[i, j] + ydrData[i, j] * ydrData[i, j] + VxxData[i, j] * VxxData[i, j] + VyyData[i, j] * VyyData[i, j] + 2.0 * VxyData[i, j] * VxyData[i, j]);
                }
            }
            return HFThetaData;
        }

        public static double GetDataMin(double[,] data)
        {
            double min = double.PositiveInfinity;
            int length0 = data.GetLength(0);
            int length1 = data.GetLength(1);
            for (int i = 0; i < length0; i++)
            {
                for (int j = 0; j < length1; j++)
                {
                    min = Math.Min(min, data[i, j]);
                }
            }
            return min;
        }

        public static double GetDataMax(double[,] data)
        {
            double max = double.NegativeInfinity;
            int length0 = data.GetLength(0);
            int length1 = data.GetLength(1);
            for (int i = 0; i < length0; i++)
            {
                for (int j = 0; j < length1; j++)
                {
                    if (data[i, j] != 1.70141e38)
                    {
                        max = Math.Max(max, data[i, j]);
                    }

                }
            }
            return max;
        }

        public static int nx, ny, ny1, ny2, nx1, nx2;
        public static double[,] expanData, boundaryX, boundaryY, cornerPoint;

        public static double[,] MinCurveExpand(double[,] data)
        {
            ny = GetExpanNum(nrows);
            nx = GetExpanNum(ncols);
            return MinCurveExpand(data, ny, nx);
        }

        public static double[,] MinCurveExpand(double[,] data, int row, int col)
        {
            ny = row;
            nx = col;
            if (ny == 0 && nx == 0)
                return null;
            if (ny == nrows && nx == ncols)
                return data;
            ny1 = (ny - nrows) / 2;
            ny2 = ny1 + nrows;
            nx1 = (nx - ncols) / 2;
            nx2 = nx1 + ncols;

            expanData = new double[ny, nx];
            boundaryX = new double[4, nx];
            boundaryY = new double[ny, 4];
            cornerPoint = new double[2, 2];

            for (int i = ny1; i < ny2; i++)
            {
                for (int j = nx1; j < nx2; j++)
                {
                    expanData[i, j] = data[i - ny1, j - nx1];
                }
            }

            var listx = new List<int>();
            var listy = new List<int>();
            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                {
                    if (expanData[i, j] == 0)
                    {
                        listy.Add(i);
                        listx.Add(j);
                    }

                }
            }

            var blackPointNum = listx.Count;
            var blankPointX = listx.ToArray();
            var blankPointY = listy.ToArray();

            CosBlankPointInitial(ny1, ny2, nx1, nx2);

            double alpha = xtick / ytick;
            double alpha2 = Math.Pow(alpha, 4);
            double alpha1 = -0.5 / (3 + 4 * alpha * alpha + 3 * alpha2);
            double alpha3 = 2.0 * alpha * alpha;
            double alpha4 = 4.0 * (1.0 + alpha * alpha);
            double alpha5 = alpha4 * alpha * alpha;

            double error = 1;
            int itera = 0;
            while (error > 1e-6 && itera < 1000)
            {
                MinCurvaExpanBoundary();
                error = 0;
                for (int n = 0; n < blackPointNum; n++)
                {
                    int i = blankPointY[n];
                    int j = blankPointX[n];

                    double temp = alpha1 * ((GetData(i + 2, j) + GetData(i - 2, j)) + alpha2 * (GetData(i, j + 2) + GetData(i, j - 2))
                        + alpha3 * (GetData(i + 1, j + 1) + GetData(i + 1, j - 1) + GetData(i - 1, j + 1) + GetData(i - 1, j - 1))
                        - alpha4 * (GetData(i + 1, j) + GetData(i - 1, j)) - alpha5 * (GetData(i, j + 1) + GetData(i, j - 1)));
                    error += Math.Abs(temp - expanData[i, j]);
                    expanData[i, j] = temp;
                }
                itera++;
                error /= blackPointNum;
            }

            return expanData;
        }

        public static double GetData(int i, int j)
        {
            if (i == -2)
            {
                return boundaryX[0, j];
            }
            else if (i == -1)
            {
                if (j == -1)
                    return cornerPoint[0, 0];
                else if (j == nx)
                    return cornerPoint[0, 1];
                else
                    return boundaryX[1, j];
            }
            else if (i == ny)
            {
                if (j == -1)
                    return cornerPoint[1, 0];
                else if (j == nx)
                    return cornerPoint[1, 1];
                else
                    return boundaryX[2, j];
            }
            else if (i == ny + 1)
            {
                return boundaryX[3, j];
            }
            else if (j == -2)
            {
                return boundaryY[i, 0];
            }
            else if (j == -1)
            {
                return boundaryY[i, 1];
            }
            else if (j == nx)
            {
                return boundaryY[i, 2];
            }
            else if (j == nx + 1)
            {
                return boundaryY[i, 3];
            }
            else
            {
                return expanData[i, j];
            }

        }

        public static void CosBlankPointInitial(int ny1, int ny2, int nx1, int nx2)
        {
            int num = nrows * 2 + (ncols - 2) * 2;
            double sum = 0.0;

            for (int i = ny1; i < ny2; i++)
            {
                sum += expanData[i, nx1] + expanData[i, nx2 - 1];
            }
            for (int j = nx1 + 1; j < nx2 - 1; j++)
            {
                sum += expanData[ny1, j] + expanData[ny2 - 1, j];
            }
            double ave = sum / num;
            for (int i = 0; i < ny; i++)
            {
                expanData[i, 0] = ave;
                expanData[i, nx - 1] = ave;
            }
            for (int j = 1; j < nx - 1; j++)
            {
                expanData[0, j] = ave;
                expanData[ny - 1, j] = ave;
            }


            var d1 = new double[ny, nx];
            var d2 = new double[ny, nx];
            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                {
                    d1[i, j] = expanData[i, j];
                    d2[i, j] = expanData[i, j];
                }
            }


            for (int i = 1; i < ny - 1; i++)
            {
                for (int j = 1; j < nx1; j++)
                {
                    d1[i, j] = (d1[ny - 1, j] - d1[0, j]) * Math.Cos((ny - 1 - i) * Math.PI / (ny - 1 - 0) / 2.0) + d1[0, j];
                }
                for (int j = nx2; j < nx - 1; j++)
                {
                    d1[i, j] = (d1[ny - 1, j] - d1[0, j]) * Math.Cos((ny - 1 - i) * Math.PI / (ny - 1 - 0) / 2.0) + d1[0, j];
                }
            }
            for (int j = nx1; j < nx2; j++)
            {
                for (int i = 1; i < ny1; i++)
                {
                    d1[i, j] = (d1[ny1, j] - d1[0, j]) * Math.Cos((ny1 - i) * Math.PI / (ny1 - 0) / 2.0) + d1[0, j];
                }
                for (int i = ny2; i < ny - 1; i++)
                {
                    d1[i, j] = (d1[ny - 1, j] - d1[ny2 - 1, j]) * Math.Cos((ny - 1 - i) * Math.PI / (ny - ny2) / 2.0) + d1[ny2 - 1, j];
                }
            }

            for (int j = 1; j < nx - 1; j++)
            {
                for (int i = 1; i < ny1; i++)
                {
                    d2[i, j] = (d2[i, nx - 1] - d2[i, 0]) * Math.Cos((nx - 1 - j) * Math.PI / (nx - 1 - 0) / 2.0) + d2[i, 0];
                }
                for (int i = ny2; i < ny - 1; i++)
                {
                    d2[i, j] = (d2[i, nx - 1] - d2[i, 0]) * Math.Cos((nx - 1 - j) * Math.PI / (nx - 1 - 0) / 2.0) + d2[i, 0];
                }
            }
            for (int i = ny1; i < ny2; i++)
            {
                for (int j = 1; j < nx1; j++)
                {
                    d2[i, j] = (d2[i, nx1] - d2[i, 0]) * Math.Cos((nx1 - j) * Math.PI / (nx1 - 0) / 2.0) + d2[i, 0];
                }
                for (int j = nx2; j < nx - 1; j++)
                {
                    d2[i, j] = (d2[i, nx - 1] - d2[i, nx2 - 1]) * Math.Cos((nx - 1 - j) * Math.PI / (nx - nx2) / 2.0) + d2[i, nx2 - 1];
                }
            }

            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                {
                    expanData[i, j] = (d1[i, j] + d2[i, j]) / 2.0;
                }
            }


        }

        public static void MinCurvaExpanBoundary()
        {
            double alpha = xtick / ytick;
            double alpha1 = alpha * alpha;
            double alpha2 = 2.0 * (1 + alpha1);
            double beta = 1.0 / alpha;
            double beta1 = beta * beta;
            double beta2 = 2.0 * (1 + beta1);

            for (int j = 0; j < nx; j++)
            {
                boundaryX[1, j] = 2.0 * expanData[0, j] - expanData[1, j];
                boundaryX[2, j] = 2.0 * expanData[ny - 1, j] - expanData[ny - 2, j];
            }

            for (int i = 0; i < ny; i++)
            {
                boundaryY[i, 1] = 2.0 * expanData[i, 0] - expanData[i, 1];
                boundaryY[i, 2] = 2.0 * expanData[i, nx - 1] - expanData[i, nx - 2];
            }

            cornerPoint[0, 0] = boundaryY[1, 1] + boundaryX[1, 1] - expanData[1, 1];
            cornerPoint[1, 0] = boundaryX[2, 1] + boundaryY[ny - 2, 1] - expanData[ny - 2, 1];
            cornerPoint[0, 1] = boundaryY[1, 2] + boundaryX[1, nx - 2] - expanData[1, nx - 2];
            cornerPoint[1, 1] = boundaryX[2, nx - 2] + boundaryY[ny - 2, 2] - expanData[ny - 2, nx - 2];

            double p;
            for (int j = 1; j < nx - 1; j++)
            {
                p = alpha1 * (expanData[1, j + 1] - boundaryX[1, j + 1] + expanData[1, j - 1] - boundaryX[1, j - 1]) - alpha2 * (expanData[1, j] - boundaryX[1, j]);
                boundaryX[0, j] = expanData[2, j] + p;
                p = alpha1 * (boundaryX[2, j + 1] - expanData[ny - 2, j + 1] + boundaryX[2, j - 1] - expanData[ny - 2, j - 1]) - alpha2 * (boundaryX[2, j] - expanData[ny - 2, j]);
                boundaryX[3, j] = expanData[ny - 3, j] - p;
            }
            //j==0
            p = alpha1 * (expanData[1, 1] - boundaryX[1, 1] + boundaryY[1, 1] - cornerPoint[0, 0]) - alpha2 * (expanData[1, 0] - boundaryX[1, 0]);
            boundaryX[0, 0] = expanData[2, 0] + p;
            p = alpha1 * (boundaryX[2, 1] - expanData[ny - 2, 1] + cornerPoint[1, 0] - boundaryY[ny - 2, 1]) - alpha2 * (boundaryX[2, 0] - expanData[ny - 2, 0]);
            boundaryX[3, 0] = expanData[ny - 3, 0] - p;
            //j==nx-1
            p = alpha1 * (boundaryY[1, 2] - cornerPoint[0, 1] + expanData[1, ny - 2] - boundaryX[1, ny - 2]) - alpha2 * (expanData[1, ny - 1] - boundaryX[1, ny - 1]);
            boundaryX[0, ny - 1] = expanData[2, ny - 1] + p;
            p = alpha1 * (cornerPoint[1, 1] - boundaryY[ny - 2, 2] + boundaryX[2, nx - 2] - expanData[ny - 2, nx - 2]) - alpha2 * (boundaryX[2, nx - 1] - expanData[ny - 2, nx - 1]);
            boundaryX[3, nx - 1] = expanData[ny - 3, nx - 1] - p;

            for (int i = 1; i < ny - 1; i++)
            {
                p = beta1 * (expanData[i + 1, 1] - boundaryY[i + 1, 1] + expanData[i - 1, 1] - boundaryY[i - 1, 1]) - beta2 * (expanData[i, 1] - boundaryY[i, 1]);
                boundaryY[i, 0] = expanData[i, 2] + p;
                p = beta1 * (boundaryY[i + 1, 2] - expanData[i + 1, nx - 2] + boundaryY[i - 1, 2] - expanData[i - 1, nx - 2]) - beta2 * (boundaryY[i, 2] - expanData[i, nx - 2]);
                boundaryY[i, 3] = expanData[i, ny - 2] - p;
            }
            //i==0
            p = beta1 * (expanData[1, 1] - boundaryY[1, 1] + boundaryX[1, 1] - cornerPoint[0, 0]) - beta2 * (expanData[0, 1] - boundaryY[0, 1]);
            boundaryY[0, 0] = expanData[0, 2] + p;
            p = beta1 * (boundaryY[1, 2] - expanData[1, nx - 2] + cornerPoint[0, 1] - boundaryX[1, nx - 2]) - beta2 * (boundaryY[0, 2] - expanData[0, nx - 2]);
            boundaryY[0, 3] = expanData[0, ny - 2] - p;
            //i=ny-1
            p = beta1 * (boundaryX[2, 1] - cornerPoint[1, 0] + expanData[ny - 2, 1] - boundaryY[ny - 2, 1]) - beta2 * (expanData[ny - 1, 1] - boundaryY[ny - 1, 1]);
            boundaryY[ny - 1, 0] = expanData[ny - 1, 2] + p;
            p = beta1 * (cornerPoint[1, 1] - boundaryX[2, nx - 2] + boundaryY[ny - 2, 2] - expanData[ny - 2, nx - 2]) - beta2 * (boundaryY[ny - 1, 2] - expanData[ny - 1, nx - 2]);
            boundaryY[ny - 1, 3] = expanData[ny - 1, ny - 2] - p;

        }

        public static int GetExpanNum(int num)
        {
            if (num <= 64)
                return 128;
            else if (num <= 128)
                return 256;
            else if (num <= 256)
                return 512;
            else if (num <= 512)
                return 1024;
            else
                return 0;
        }

        public static double[,] ExpandRestore(double[,] expanData)
        {
            if (ny == nrows && nx == ncols)
            {
                return expanData;
            }
            var data = new double[nrows, ncols];
            for (int i = ny1; i < ny2; i++)
            {
                for (int j = nx1; j < nx2; j++)
                {
                    data[i - ny1, j - nx1] = expanData[i, j];
                }
            }
            return data;
        }

        public static double[,] GuassNoise(double[,] data)
        {
            var random = new Random(9999);
            double r1, r2, noise;
            int nrows = data.GetLength(0);
            int ncols = data.GetLength(1);
            var result = new double[nrows, ncols];

            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    r1 = random.NextDouble();
                    r2 = random.NextDouble();
                    noise = Math.Sqrt((-2) * Math.Log(r2)) * Math.Sin(2 * Math.PI * r1);
                    result[i, j] = data[i, j] * (1 + noise * 0.05);
                }
            }
            return result;
        }
    }
}
