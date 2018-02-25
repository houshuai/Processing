using System;
using System.Windows.Media;

namespace Processing
{
    public class Colormap
    {
        public int ColormapLength, Ydivisions;
        public byte AlphaValue;
        public ColormapBrushEnum ColormapBrushType;
        public Color[] colors;

        public Colormap()
        {
            ColormapLength = 64;
            AlphaValue = 255;
            Ydivisions = 64;
            ColormapBrushType = ColormapBrushEnum.Jet;
        }

        public enum ColormapBrushEnum
        {
            Spring = 0,
            Summer = 1,
            Autumn = 2,
            Winter = 3,
            Gray = 4,
            Jet = 5,
            Hot = 6,
            Cool = 7
        }

        public void SetColors()
        {
            byte[,] cmap = new byte[ColormapLength, 4];
            double[] array = new double[ColormapLength];

            switch (ColormapBrushType)
            {
                case ColormapBrushEnum.Spring:
                    for (int i = 0; i < ColormapLength; i++)
                    {
                        array[i] = 1.0 * i / (ColormapLength - 1);
                        cmap[i, 0] = AlphaValue;
                        cmap[i, 1] = 255;
                        cmap[i, 2] = (byte)(255 * array[i]);
                        cmap[i, 3] = (byte)(255 - cmap[i, 2]);
                    }
                    break;
                case ColormapBrushEnum.Summer:
                    for (int i = 0; i < ColormapLength; i++)
                    {
                        array[i] = 1.0 * i / (ColormapLength - 1);
                        cmap[i, 0] = AlphaValue;
                        cmap[i, 1] = (byte)(255 * array[i]);
                        cmap[i, 2] = (byte)(255 * 0.5 * (1 + array[i]));
                        cmap[i, 3] = (byte)(255 * 0.4);
                    }
                    break;
                case ColormapBrushEnum.Autumn:
                    for (int i = 0; i < ColormapLength; i++)
                    {
                        array[i] = 1.0 * i / (ColormapLength - 1);
                        cmap[i, 0] = AlphaValue;
                        cmap[i, 1] = 255;
                        cmap[i, 2] = (byte)(255 * array[i]);
                        cmap[i, 3] = 0;
                    }
                    break;
                case ColormapBrushEnum.Winter:
                    for (int i = 0; i < ColormapLength; i++)
                    {
                        array[i] = 1.0 * i / (ColormapLength - 1);
                        cmap[i, 0] = AlphaValue;
                        cmap[i, 1] = 0;
                        cmap[i, 2] = (byte)(255 * array[i]);
                        cmap[i, 3] = (byte)(255 * (1.0 - 0.5 * cmap[i, 2]));
                    }
                    break;
                case ColormapBrushEnum.Gray:
                    for (int i = 0; i < ColormapLength; i++)
                    {
                        array[i] = 1.0 * i / (ColormapLength - 1);
                        cmap[i, 0] = AlphaValue;
                        cmap[i, 1] = (byte)(255 * array[i]);
                        cmap[i, 2] = (byte)(255 * array[i]);
                        cmap[i, 3] = (byte)(255 * array[i]);
                    }
                    break;
                case ColormapBrushEnum.Jet:
                    int n = (int)Math.Ceiling(ColormapLength / 4.0);
                    double[,] cMatrix = new double[ColormapLength, 3];
                    int nMod = 0;
                    double[] array1 = new double[3 * n - 1];
                    int[] red = new int[array1.Length];
                    int[] green = new int[array1.Length];
                    int[] blue = new int[array1.Length];

                    if (ColormapLength % 4 == 1)
                        nMod = 1;

                    for (int i = 0; i < array1.Length; i++)
                    {
                        if (i < n)
                        {
                            array1[i] = (i + 1.0) / n;
                        }
                        else if (i >= n && i < 2 * n - 1)
                        {
                            array1[i] = 1.0;
                        }
                        else if (i >= 2 * n - 1)
                        {
                            array1[i] = (3.0 * n - 1.0 - i) / n;
                        }
                        green[i] = (int)Math.Ceiling(n / 2.0) - nMod + i;
                        red[i] = green[i] + n;
                        blue[i] = green[i] - n;
                    }

                    int nb = 0;
                    for (int i = 0; i < blue.Length; i++)
                    {
                        if (blue[i] > 0)
                            nb++;
                    }
                    for (int i = 0; i < ColormapLength; i++)
                    {
                        for (int j = 0; j < red.Length; j++)
                        {
                            if (i == red[j] && red[j] < ColormapLength)
                                cMatrix[i, 0] = array1[i - red[0]];
                        }
                        for (int j = 0; j < green.Length; j++)
                        {
                            if (i == green[j] && green[j] < ColormapLength)
                                cMatrix[i, 1] = array1[i - green[0]];
                        }
                        for (int j = 0; j < blue.Length; j++)
                        {
                            if (i == blue[j] && blue[j] >= 0)
                                cMatrix[i, 2] = array1[array1.Length - 1 - nb + i];
                        }
                    }
                    for (int i = 0; i < ColormapLength; i++)
                    {
                        cmap[i, 0] = AlphaValue;
                        for (int j = 0; j < 3; j++)
                        {
                            cmap[i, j + 1] = (byte)(cMatrix[i, j] * 255);
                        }
                    }
                    break;
                case ColormapBrushEnum.Hot:
                    int n1 = 3 * ColormapLength / 8;
                    double[] red1 = new double[ColormapLength];
                    double[] green1 = new double[ColormapLength];
                    double[] blue1 = new double[ColormapLength];
                    for (int i = 0; i < ColormapLength; i++)
                    {
                        if (i < n1)
                        {
                            red1[i] = 1.0 * (i + 1.0) / n1;
                        }
                        else
                        {
                            red1[i] = 1.0;
                        }
                        if (i < n1)
                        {
                            green1[i] = 0.0;
                        }
                        else if (i >= n1 && i < 2 * n1)
                        {
                            green1[i] = 1.0 * (i + 1 - n1) / n1;
                        }
                        else
                        {
                            green1[i] = 1.0;
                        }
                        if (i < 2 * n1)
                        {
                            blue1[i] = 0.0;
                        }
                        else
                        {
                            blue1[i] = 1.0 * (i + 1 - 2 * n1) / (ColormapLength - 2 * n1);
                        }

                        cmap[i, 0] = AlphaValue;
                        cmap[i, 1] = (byte)(255 * red1[i]);
                        cmap[i, 2] = (byte)(255 * green1[i]);
                        cmap[1, 3] = (byte)(255 * blue1[i]);
                    }
                    break;
                case ColormapBrushEnum.Cool:
                    for (int i = 0; i < ColormapLength; i++)
                    {
                        array[i] = 1.0 * i / (ColormapLength - 1);
                        cmap[i, 0] = AlphaValue;
                        cmap[i, 1] = (byte)(255 * array[i]);
                        cmap[i, 2] = (byte)(255 * (1 - array[i]));
                        cmap[i, 3] = 255;
                    }
                    break;
                default:
                    break;
            }

            colors = new Color[Ydivisions];
            for (int i = 0; i < Ydivisions; i++)
            {
                int colorIndex = (ColormapLength - 1) * i / (Ydivisions - 1);
                colors[i] = Color.FromArgb(
                    cmap[colorIndex, 0],
                    cmap[colorIndex, 1],
                    cmap[colorIndex, 2],
                    cmap[colorIndex, 3]);
            }
        }
        
        public Color GetColor(double y,double ymin,double ymax)
        {
            double dy = (ymax - ymin) / (Ydivisions - 1);
            for (int i = 0; i < Ydivisions; i++)
            {
                double y1 = ymin + i * dy;
                if (y >= y1 && y < y1 + dy)
                {
                    return colors[i];
                }
            }
            return Colors.Transparent;
        }
    }
}
