using System;
using MathNet.Numerics.LinearAlgebra;
using System.Collections.Generic;
using System.Windows;

namespace Processing
{
    public static class Separation
    {
        private static MatrixBuilder<double> Matrix = Matrix<double>.Build;
        private static VectorBuilder<double> Vector = Vector<double>.Build;

        public static double[,] TrendAnalysis(double[,] z, int pow)
        {
            int row = z.GetLength(0);
            int col = z.GetLength(1);
            var n = (pow + 1) * (pow + 2) / 2;
            //var m = row * col;

            var result = new double[row, col];

            var zv = new List<double>();
            var corev = new List<double[]>();
            int currXPow, currYPow, powIndex;
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    if (z[i, j] != 1.70141e38)
                    {
                        zv.Add(z[i, j]);
                        var co = new double[n];
                        powIndex = currXPow = currYPow = 0;
                        for (int k = 0; k < n; k++)
                        {
                            co[k] = Math.Pow(i + 1, currXPow) * Math.Pow(j + 1, currYPow);
                            if (currYPow == powIndex)
                            {
                                powIndex++;
                                currXPow = powIndex;
                                currYPow = 0;
                            }
                            else
                            {
                                currXPow--;
                                currYPow++;
                            }
                        }
                        corev.Add(co);
                    }
                }
            }

            int count = corev.Count;
            var coreArray = new double[count, n];
            for (int i = 0; i < count; i++)
            {
                var co = corev[i];
                for (int j = 0; j < n; j++)
                {
                    coreArray[i, j] = co[j];
                }
            }

            var zx = Vector.Dense(zv.ToArray());
            var core = Matrix.DenseOfArray(coreArray);
            var coreT = core.Transpose();
            var temp = coreT * core;
            var a = temp.Inverse() * coreT * zx;
            var rv = (Matrix.DenseOfArray(GetCoreMatrix(row * col, n, col)) * a).AsArray();

            int index = 0;
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    result[i, j] = rv[index];
                    index++;
                }
            }

            return result;
        }

        private static double[,] GetCoreMatrix(int m, int n, int col)
        {
            int x = 0, y = 0;
            var core = new double[m, n];

            int currXPow, currYPow, powIndex;
            for (int i = 0; i < m; i++)
            {
                if (i % col == 0)
                {
                    x++;
                    y = 1;
                }
                else
                {
                    y++;
                }
                powIndex = currXPow = currYPow = 0;
                for (int j = 0; j < n; j++)
                {
                    core[i, j] = Math.Pow(x, currXPow) * Math.Pow(y, currYPow);
                    if (currYPow == powIndex)
                    {
                        powIndex++;
                        currXPow = powIndex;
                        currYPow = 0;
                    }
                    else
                    {
                        currXPow--;
                        currYPow++;
                    }
                }
            }

            return core;
        }

        public static double[,] FilterTrendAnalysis(double[,] z, double[,] filted, int pow)
        {
            return TrendAnalysis(filted, pow);
        }

        public static double[,] Filter(double[,] z, int expendMethod, double filterProject, int length)
        {
            int row = z.GetLength(0);
            int col = z.GetLength(1);

            int expendRow = EdgeDetection.GetExpanNum(row);
            int expendCol = EdgeDetection.GetExpanNum(col);
            int row1 = (expendRow - row) / 2;
            int row2 = expendRow - row1 - row;
            int col1 = (expendCol - col) / 2;
            int col2 = expendCol - col1 - col;

            var expend = Expend(z, expendMethod, row1, row2, col1, col2);
            var frqData = FFT.fft_2D(new Complex2D(expend));
            var zdrfData = ExpandRestore(FFT.idft_2D(EdgeDetection.Zdrf(frqData)).GetRe(), row, col, row1, col1);

            return Filter(z, zdrfData, filterProject, length);
        }

        public static double[,] Filter(double[,] z, double[,] zdr, double filterProject, int length)
        {
            int row = z.GetLength(0);
            int col = z.GetLength(1);
            var project = GetProject(zdr);

            var result = new double[row, col];
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    result[i, j] = z[i, j];
                }
            }
            var isReached = new bool[row, col];
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    if (project[i, j] < filterProject)
                    {
                        isReached[i, j] = true;
                        //project[i, j] = 1.70141e38;
                    }
                }
            }
            for (int i = 0; i < row; i++)
            {
                isReached[i, 0] = true;
                isReached[i, col - 1] = true;
            }
            for (int j = 0; j < col; j++)
            {
                isReached[0, j] = true;
                isReached[row - 1, j] = true;
            }

            DeleteLocalRegion(result, isReached, length, project, filterProject);

            return result;
        }

        public static double[,] FilterProject(double[,] z, int expendMethod, double filterProject)
        {
            int row = z.GetLength(0);
            int col = z.GetLength(1);

            int expendRow = EdgeDetection.GetExpanNum(row);
            int expendCol = EdgeDetection.GetExpanNum(col);
            int row1 = (expendRow - row) / 2;
            int row2 = expendRow - row1 - row;
            int col1 = (expendCol - col) / 2;
            int col2 = expendCol - col1 - col;

            var expend = Expend(z, expendMethod, row1, row2, col1, col2);
            var frqData = FFT.fft_2D(new Complex2D(expend));
            var zdrfData = ExpandRestore(FFT.idft_2D(EdgeDetection.Zdrf(frqData)).GetRe(), row, col, row1, col1);

            var project = GetProject(zdrfData);

            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    if (project[i, j] < filterProject)
                    {
                        project[i, j] = 1.70141e38;
                    }
                }
            }

            return project;
        }

        private static double[,] GetProject(double[,] zdr)
        {
            int row = zdr.GetLength(0);
            int col = zdr.GetLength(1);

            var normal = Normal(zdr);
            var up = new Vector3D(0, 0, 1);
            var project = new double[row, col];
            var array = new double[row * col];
            int index = 0;
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    project[i, j] = Vector3D.Dot(normal[i, j], up);
                    array[index++] = project[i, j];
                }
            }

            int rangeNum = 10;
            double rangeTick = 1.0 / rangeNum;

            QuickSort(array, 0, array.Length - 1);

            int countOfOneRange = array.Length / rangeNum;
            var range = new double[rangeNum];
            for (int i = 0; i < rangeNum; i++)
            {
                range[i] = array[i * countOfOneRange];
            }


            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    var value = project[i, j];
                    for (int k = rangeNum - 1; k >= 0; k--)
                    {
                        if (value > range[k])
                        {
                            double oldMin = range[k];
                            double oldMax;
                            if (k == rangeNum - 1)
                            {
                                oldMax = array[array.Length - 1];
                            }
                            else
                            {
                                oldMax = range[k + 1];
                            }
                            project[i, j] = (project[i, j] - oldMin) / (oldMax - oldMin) * rangeTick + rangeTick * k;
                            break;
                        }
                    }
                }
            }

            return project;
        }

        private static void QuickSort(double[] array, int low, int high)
        {
            if (low < high)
            {
                int left = low;
                int right = high;
                double temp = array[low];

                while (low < high)
                {
                    while (low < high && array[high] >= temp)
                    {
                        high--;
                    }
                    array[low] = array[high];
                    while (low < high && array[low] <= temp)
                    {
                        low++;
                    }
                    array[high] = array[low];
                }

                array[low] = temp;

                QuickSort(array, left, low - 1);
                QuickSort(array, low + 1, right);
            }
        }

        private static void DeleteLocalRegion(double[,] data, bool[,] isSearched, int length, double[,] project, double filterProject)
        {
            int row = data.GetLength(0);
            int col = data.GetLength(1);

            var localSearched = new bool[row, col];
            for (int i = 0; i < row; i++)
            {
                localSearched[i, 0] = true;
                localSearched[i, col - 1] = true;
            }
            for (int j = 0; j < col; j++)
            {
                localSearched[0, j] = true;
                localSearched[row - 1, j] = true;
            }

            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    if (!isSearched[i, j])
                    {
                        //获得去除小于filterProject包围的局部区域
                        RectInt rect;
                        var pointList = SmallRegion(i, j, isSearched, length, out rect);

                        if (pointList != null)
                        {
                            //获取局部区域
                            pointList = LocalRegion(pointList, rect, localSearched, project, filterProject);
                        }

                        if (pointList != null)
                        {
                            //删除点
                            foreach (var point in pointList)
                            {
                                data[point.X, point.Y] = 1.70141e38;
                            }
                        }
                    }

                }
            }
        }

        private static List<PointInt> SmallRegion(int x, int y, bool[,] isSearched, int length, out RectInt rect)
        {
            rect = new RectInt(x, y, x, y);
            var pointList = new List<PointInt>();
            var queue = new Queue<PointInt>();
            queue.Enqueue(new PointInt(x, y));
            isSearched[x, y] = true;

            while (queue.Count > 0)
            {
                var point = queue.Dequeue();

                rect.X1 = rect.X1 < point.X ? rect.X1 : point.X;
                rect.Y1 = rect.Y1 < point.Y ? rect.Y1 : point.Y;
                rect.X2 = rect.X2 > point.X ? rect.X2 : point.X;
                rect.Y2 = rect.Y2 > point.Y ? rect.Y2 : point.Y;
                pointList.Add(point);

                if (!isSearched[point.X - 1, point.Y])
                {
                    queue.Enqueue(new PointInt(point.X - 1, point.Y));
                    isSearched[point.X - 1, point.Y] = true;
                }
                if (!isSearched[point.X + 1, point.Y])
                {
                    queue.Enqueue(new PointInt(point.X + 1, point.Y));
                    isSearched[point.X + 1, point.Y] = true;
                }
                if (!isSearched[point.X, point.Y - 1])
                {
                    queue.Enqueue(new PointInt(point.X, point.Y - 1));
                    isSearched[point.X, point.Y - 1] = true;
                }
                if (!isSearched[point.X, point.Y + 1])
                {
                    queue.Enqueue(new PointInt(point.X, point.Y + 1));
                    isSearched[point.X, point.Y + 1] = true;
                }

            }

            if (rect.MaxSide < length)
            {
                return pointList;
            }

            return null;
        }

        private static List<PointInt> LocalRegion(List<PointInt> smallList, RectInt smallRect, bool[,] isSearched, double[,] project, double filterProject)
        {
            var pointList = new List<PointInt>();
            var queue = new Queue<PointInt>();
            foreach (var point in smallList)
            {
                queue.Enqueue(point);
                isSearched[point.X, point.Y] = true;
            }
            var rect = new RectInt(smallRect.X1 - smallRect.XLen / 2, smallRect.Y1 - smallRect.YLen / 2,
                smallRect.X2 + smallRect.XLen / 2, smallRect.Y2 + smallRect.YLen / 2);

            while (queue.Count > 0)
            {
                var point = queue.Dequeue();
                var filter = project[point.X, point.Y];
                pointList.Add(point);

                var newPoint = new PointInt(point.X - 1, point.Y);
                Test(newPoint, isSearched, project, queue, filterProject);
                newPoint = new PointInt(point.X + 1, point.Y);
                Test(newPoint, isSearched, project, queue, filterProject);
                newPoint = new PointInt(point.X, point.Y - 1);
                Test(newPoint, isSearched, project, queue, filterProject);
                newPoint = new PointInt(point.X, point.Y + 1);
                Test(newPoint, isSearched, project, queue, filterProject);

            }

            return pointList;
        }

        private static void Test(PointInt newPoint, bool[,] isSearched, double[,] project, Queue<PointInt> queue, double filter)
        {
            var newProject = project[newPoint.X, newPoint.Y];
            if (!isSearched[newPoint.X, newPoint.Y] && newProject < filter)
            {
                queue.Enqueue(newPoint);
                isSearched[newPoint.X, newPoint.Y] = true;
            }
        }

        private static void Delete(double[,] project, int x, int y, double[,] result, double filterProject)
        {
            result[x, y] = 1.70141e38;
            if (x > 1 && result[x - 1, y] != 1.70141e38 && project[x - 1, y] < filterProject)
            {
                Delete(project, x - 1, y, result, filterProject);
            }

            if (y > 1 && result[x, y - 1] != 1.70141e38 && project[x, y - 1] < filterProject)
            {
                Delete(project, x, y - 1, result, filterProject);
            }

            if (x < project.GetLength(0) - 2 && result[x + 1, y] != 1.70141e38 && project[x + 1, y] < filterProject)
            {
                Delete(project, x + 1, y, result, filterProject);
            }

            if (y < project.GetLength(1) - 2 && result[x, y + 1] != 1.70141e38 && project[x, y + 1] < filterProject)
            {
                Delete(project, x, y + 1, result, filterProject);
            }
        }

        private static Vector3D[,] Normal(double[,] z)
        {
            int row = z.GetLength(0);
            int col = z.GetLength(1);

            var normal = new Vector3D[row, col];

            for (int i = 1; i < row - 1; i++)
            {
                for (int j = 1; j < col - 1; j++)
                {
                    var v1 = new Vector3D(1, 0, z[i + 1, j] - z[i, j]);
                    var v2 = new Vector3D(0, 1, z[i, j + 1] - z[i, j]);
                    var v3 = new Vector3D(-1, 0, z[i - 1, j] - z[i, j]);
                    var v4 = new Vector3D(0, -1, z[i, j - 1] - z[i, j]);
                    var n1 = Vector3D.Cross(v1, v2);
                    var n2 = Vector3D.Cross(v2, v3);
                    var n3 = Vector3D.Cross(v3, v4);
                    var n4 = Vector3D.Cross(v4, v1);
                    var n = n1 + n2 + n3 + n4;
                    normal[i, j] = n.Normalize();
                }
            }

            for (int i = 1; i < row - 1; i++)
            {
                var v1 = new Vector3D(1, 0, z[i + 1, 0] - z[i, 0]);
                var v2 = new Vector3D(0, 1, z[i, 1] - z[i, 0]);
                var v3 = new Vector3D(-1, 0, z[i - 1, 0] - z[i, 0]);
                var n1 = Vector3D.Cross(v1, v2);
                var n2 = Vector3D.Cross(v2, v3);
                var n = n1 + n2;
                normal[i, 0] = n.Normalize();

                var v11 = new Vector3D(1, 0, z[i + 1, col - 1] - z[i, col - 1]);
                var v31 = new Vector3D(-1, 0, z[i - 1, col - 1] - z[i, col - 1]);
                var v4 = new Vector3D(0, -1, z[i, col - 2] - z[i, col - 1]);
                var n3 = Vector3D.Cross(v31, v4);
                var n4 = Vector3D.Cross(v4, v11);
                n = n3 + n4;
                normal[i, col - 1] = n.Normalize();
            }

            for (int j = 1; j < col - 1; j++)
            {
                var v1 = new Vector3D(1, 0, z[1, j] - z[0, j]);
                var v2 = new Vector3D(0, 1, z[0, j + 1] - z[0, j]);
                var v4 = new Vector3D(0, -1, z[0, j - 1] - z[0, j]);
                var n1 = Vector3D.Cross(v1, v2);
                var n4 = Vector3D.Cross(v4, v1);
                var n = n1 + n4;
                normal[0, j] = n.Normalize();

                var v21 = new Vector3D(0, 1, z[row - 1, j + 1] - z[row - 1, j]);
                var v3 = new Vector3D(-1, 0, z[row - 2, j] - z[row - 1, j]);
                var v41 = new Vector3D(0, -1, z[row - 1, j - 1] - z[row - 1, j]);
                var n2 = Vector3D.Cross(v21, v3);
                var n3 = Vector3D.Cross(v3, v41);
                n = n2 + n3;
                normal[row - 1, j] = n.Normalize();
            }

            var n0 = Vector3D.Cross(new Vector3D(1, 0, z[1, 0] - z[0, 0]), new Vector3D(0, 1, z[0, 1] - z[0, 0]));
            normal[0, 0] = n0.Normalize();

            n0 = Vector3D.Cross(new Vector3D(0, 1, z[row - 1, 1] - z[row - 1, 0]), new Vector3D(-1, 0, z[row - 2, 0] - z[row - 1, 0]));
            normal[row - 1, 0] = n0.Normalize();

            n0 = Vector3D.Cross(new Vector3D(-1, 0, z[row - 2, col - 1] - z[row - 1, col - 1]), new Vector3D(0, -1, z[row - 1, col - 2] - z[row - 1, col - 1]));
            normal[row - 1, col - 1] = n0.Normalize();

            n0 = Vector3D.Cross(new Vector3D(0, -1, z[0, col - 2] - z[0, col - 1]), new Vector3D(1, 0, z[1, col - 1] - z[0, col - 1]));
            normal[0, col - 1] = n0.Normalize();

            return normal;
        }

        public static double[,] InterpolateCut(double[,] data, int expendMethod, int r, int times)
        {
            int row = data.GetLength(0);
            int col = data.GetLength(1);

            var region = data;

            int time = 0;
            double delta = 0;

            do
            {
                var origin = region;
                double[,] g = Expend(data, expendMethod, r, r, r, r);
                region = new double[row, col];
                delta = 0;
                for (int i = r; i < row + r; i++)
                {
                    for (int j = r; j < col + r; j++)
                    {
                        var tempg = g[i, j];
                        var tempgx_r = g[i - r, j];
                        var tempgy_r = g[i, j - r];
                        var tempgyr = g[i, j + r];
                        var tempgxr = g[i + r, j];
                        var deltaxa = tempg - (tempgx_r + tempgxr) / 2;
                        var deltaya = tempg - (tempgy_r + tempgyr) / 2;
                        var deltax = tempgxr - tempgx_r;
                        var deltay = tempgyr - tempgy_r;

                        double a, b, c;
                        if (deltax == 0 && deltaxa == 0)
                        {
                            b = 1;
                        }
                        else
                        {
                            b = deltax * deltax / (deltaxa * deltaxa + deltax * deltax);
                        }

                        if (deltay == 0 && deltaya == 0)
                        {
                            c = 1;
                        }
                        else
                        {
                            c = deltay * deltay / (deltaya * deltaya + deltay * deltay);
                        }


                        a = (b + c) / 2.0;
                        var value = (1.0 - a) * (tempgx_r + tempgy_r + tempgyr + tempgxr) / 4.0 + tempg * a;
                        region[i - r, j - r] = value;
                        delta += Math.Abs(value - tempg);
                    }
                }
                time++;
                delta /= row * col;

            } while (delta > 0.00001 && time < times);

            return region;
        }

        public static double[,] InterpolateCut(double[,] data, double[,] local, int expendMethod, int rMin, int rMax, int tick, int times, out int radius)
        {
            int row = data.GetLength(0);
            int col = data.GetLength(1);
            int count = (rMax - rMin) / tick;
            int r = rMin;
            radius = r;
            double delta = 0;
            double minDelta = double.MaxValue;
            double[,] result = new double[0, 0];
            int n2 = row * col * row * col;

            for (int i = 0; i < count; i++)
            {
                var region = InterpolateCut(data, expendMethod, r, times);
                delta = 0;
                for (int x = 0; x < row; x++)
                {
                    for (int y = 0; y < col; y++)
                    {
                        var s = data[x, y] - region[x, y] - local[x, y];
                        s = s * s / n2;
                        delta += s;
                    }
                }
                if (delta < minDelta)
                {
                    minDelta = delta;
                    radius = r;
                    result = region;
                }
                r += tick;
            }

            return result;
        }

        public static double[,] MoveAverage(double[,] data, int expendMethod, int r)
        {
            int row = data.GetLength(0);
            int col = data.GetLength(1);

            double[,] origin = Expend(data, expendMethod, r, r, r, r);

            var result = new double[row, col];

            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    result[i, j] = GetAverage(origin, i, j, r);
                }
            }

            return result;
        }

        public static double[,] FilterMoveAverage(double[,] data, double[,] filted, int expendMethod, int r)
        {
            var fill = AverageFill(filted);
            return MoveAverage(fill, expendMethod, r);
        }

        public static double[,] MoveAverage(double[,] data, double[,] local, int expendMethod, int rMin, int rMax, int tick, out int radius)
        {
            int row = data.GetLength(0);
            int col = data.GetLength(1);
            int count = (rMax - rMin) / tick;
            int r = rMin;
            radius = r;
            double delta = 0;
            double minDelta = double.MaxValue;
            double[,] result = new double[0, 0];
            int n2 = row * col * row * col;

            for (int i = 0; i < count; i++)
            {
                var region = MoveAverage(data, expendMethod, r);
                delta = 0;
                for (int x = 0; x < row; x++)
                {
                    for (int y = 0; y < col; y++)
                    {
                        var s = data[x, y] - region[x, y] - local[x, y];
                        s = s * s / n2;
                        delta += s;
                    }
                }
                if (delta < minDelta)
                {
                    minDelta = delta;
                    radius = r;
                    result = region;
                }
                r += tick;
            }

            return result;
        }

        private static double GetAverage(double[,] data, int xStart, int yStart, int r)
        {
            double average = 0;
            int size = 2 * r + 1;

            for (int i = xStart; i < xStart + size; i++)
            {
                for (int j = yStart; j < yStart + size; j++)
                {
                    average += data[i, j];
                }
            }

            return average / size / size;
        }

        private static double[,] Expend(double[,] data, int expendMethod, int row1, int row2, int col1, int col2)
        {
            double[,] result;
            if (expendMethod == 0)
            {
                result = GradientExpend(data, row1, row2, col1, col2);
            }
            else if (expendMethod == 1)
            {
                result = AverageExpend(data, row1, row2, col1, col2);
            }
            else
            {
                throw new InvalidOperationException("expendMethod is not exist");
            }
            return result;
        }

        private static double[,] GradientExpend(double[,] data, int row1, int row2, int col1, int col2)
        {
            int row = data.GetLength(0);
            int col = data.GetLength(1);
            var result = new double[row + row1 + row2, col + col1 + col2];

            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    result[i + row1, j + col1] = data[i, j];
                }
            }

            //left
            for (int i = 0; i < row; i++)
            {
                double gradient = 0;
                for (int j = 0; j < col1; j++)
                {
                    gradient += data[i, j] - data[i, j + 1];
                }
                gradient /= col1;

                for (int j = 0; j < col1; j++)
                {
                    result[i + row1, j] = data[i, 0] + gradient * (col1 - j);
                }
            }

            //right
            for (int i = 0; i < row; i++)
            {
                double gradient = 0;
                for (int j = col - col2; j < col; j++)
                {
                    gradient += data[i, j] - data[i, j - 1];
                }
                gradient /= col2;

                for (int j = col + col1; j < col + col1 + col2; j++)
                {
                    result[i + row1, j] = data[i, col - 1] + gradient * (j - col - col1 + 1);
                }
            }

            //top
            for (int j = 0; j < col + col1 + col2; j++)
            {
                double gradient = 0;
                for (int i = row1; i < 2 * row1; i++)
                {
                    gradient += result[i, j] - result[i + 1, j];
                }
                gradient /= row1;

                for (int i = 0; i < row1; i++)
                {
                    result[i, j] = result[row1, j] + gradient * (row1 - i);
                }
            }

            //bottom
            for (int j = 0; j < col + col1 + col2; j++)
            {
                double gradient = 0;
                for (int i = row + row1 - row2; i < row + row1; i++)
                {
                    gradient += result[i, j] - result[i - 1, j];
                }
                gradient /= row2;

                for (int i = row + row1; i < row + row1 + row2; i++)
                {
                    result[i, j] = result[row + row1 - 1, j] + gradient * (i - row - row1 + 1);
                }
            }


            return result;
        }

        private static double[,] AverageExpend(double[,] data, int row1, int row2, int col1, int col2)
        {
            int row = data.GetLength(0);
            int col = data.GetLength(1);

            var expendRow = row + row1 + row2;
            var expendCol = col + col1 + col2;
            var result = new double[expendRow, expendCol];

            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    result[i + row1, j + col1] = data[i, j];
                }
            }

            double boundAverage = 0;
            for (int i = 0; i < row; i++)
            {
                boundAverage += data[i, 0] + data[i, col - 1];
            }
            for (int j = 1; j < col - 1; j++)
            {
                boundAverage += data[0, j] + data[row - 1, j];
            }
            boundAverage /= 2 * row + 2 * (col - 2);

            for (int i = 0; i < expendRow; i++)
            {
                result[i, 0] = result[i, expendCol - 1] = boundAverage;
            }
            for (int j = 0; j < expendCol; j++)
            {
                result[0, j] = result[expendRow - 1, j] = boundAverage;
            }

            double delta = 0;
            do
            {
                delta = 0;
                for (int i = 1; i < expendRow - 1; i++)
                {
                    for (int j = 1; j < expendCol - 1; j++)
                    {
                        if ((i < row1 || i >= row1 + row) || (j < col1 || j >= col1 + col))
                        {
                            var ave = (result[i - 1, j] + result[i + 1, j] + result[i, j - 1] + result[i, j + 1]) / 4;
                            delta += Math.Abs(result[i, j] - ave);
                            result[i, j] = ave;
                        }
                    }
                }
                delta /= expendRow * expendCol - row * col;
            } while (delta > 1.0e-6);

            return result;
        }

        public static double[,] AverageFill(double[,] data)
        {
            int row = data.GetLength(0);
            int col = data.GetLength(1);
            var result = new double[row, col];
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    result[i, j] = data[i, j];
                }
            }

            double boundAverage = 0;
            int count = 0;
            for (int i = 0; i < row; i++)
            {
                if (data[i, 0] != 1.70141e38)
                {
                    boundAverage += data[i, 0];
                    count++;
                }
                if (data[i, col - 1] != 1.70141e38)
                {
                    boundAverage += data[i, col - 1];
                    count++;
                }
            }
            for (int j = 1; j < col - 1; j++)
            {
                if (data[0, j] != 1.70141e38)
                {
                    boundAverage += data[0, j];
                    count++;
                }
                if (data[row - 1, j] != 1.70141e38)
                {
                    boundAverage += data[row - 1, j];
                    count++;
                }
            }
            boundAverage /= count;

            for (int i = 0; i < row; i++)
            {
                if (result[i, 0] == 1.70141e38)
                {
                    result[i, 0] = boundAverage;
                }
                if (result[i, col - 1] == 1.70141e38)
                {
                    result[i, col - 1] = boundAverage;
                }
            }
            for (int j = 1; j < col - 1; j++)
            {
                if (result[0, j] == 1.70141e38)
                {
                    result[0, j] = boundAverage;
                }
                if (result[row - 1, j] == 1.70141e38)
                {
                    result[row - 1, j] = boundAverage;
                }
            }

            var blankList = new List<PointInt>();
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    if (result[i, j] == 1.70141e38)
                    {
                        blankList.Add(new PointInt(i, j));
                        result[i, j] = 0;
                    }
                }
            }

            double delta = 0;
            do
            {
                foreach (var blank in blankList)
                {
                    var fill = (result[blank.X - 1, blank.Y] + result[blank.X + 1, blank.Y] +
                        result[blank.X, blank.Y - 1] + result[blank.X, blank.Y + 1]) / 4;
                    delta += Math.Abs(result[blank.X, blank.Y] - fill);
                    result[blank.X, blank.Y] = fill;
                }
                delta /= blankList.Count;
            } while (delta > 1.0e-6);

            return result;
        }

        public static double[,] ExpandRestore(double[,] expanData, int row, int col, int row1, int col1)
        {
            var data = new double[row, col];
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < col; j++)
                {
                    data[i, j] = expanData[i + row1, j + col1];
                }
            }
            return data;
        }
    }

    public struct Vector3D
    {
        double X, Y, Z;

        public double Length
        {
            get
            {
                return Math.Sqrt(X * X + Y * Y + Z * Z);
            }
        }

        public Vector3D(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public static Vector3D Cross(Vector3D v1, Vector3D v2)
        {
            return new Vector3D(v1.Y * v2.Z - v1.Z * v2.Y, v1.Z * v2.X - v1.X * v2.Z, v1.X * v2.Y - v1.Y * v2.X);
        }

        public static double Dot(Vector3D v1, Vector3D v2)
        {
            return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
        }

        public static Vector3D operator +(Vector3D v1, Vector3D v2)
        {
            return new Vector3D(v1.X + v2.X, v1.Y + v2.Y, v1.Z + v2.Z);
        }

        public Vector3D Normalize()
        {
            var length = Math.Sqrt(X * X + Y * Y + Z * Z);
            return new Vector3D(X / length, Y / length, Z / length);
        }

        public Vector3D Scale(double s)
        {
            return new Vector3D(X * s, Y * s, Z * s);
        }
    }

    public struct PointInt
    {
        public int X, Y;

        public PointInt(int x, int y)
        {
            X = x;
            Y = y;
        }
    }

    public struct RectInt
    {
        public int X1, Y1, X2, Y2;
        public int XLen { get { return X2 - X1; } }
        public int YLen { get { return Y2 - Y1; } }
        public int MaxSide { get { return Math.Max(X2 - X1, Y2 - Y1); } }

        public RectInt(int x1, int y1, int x2, int y2)
        {
            X1 = x1;
            Y1 = y1;
            X2 = x2;
            Y2 = y2;
        }

        public bool Contains(PointInt point)
        {
            if (point.X > X1 && point.Y > Y1 && point.X < X2 && point.Y < Y2)
            {
                return true;
            }

            return false;
        }
    }

}
