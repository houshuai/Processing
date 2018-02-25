using System;
using System.Globalization;
using System.Windows;
using System.Windows.Data;

namespace Processing
{
    class STVConverter : IValueConverter
    {
        public object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            int selectedIndex = (int)value;
            int para = System.Convert.ToInt32(parameter);

            if (selectedIndex == para)
            {
                return Visibility.Visible;
            }
            else
            {
                return Visibility.Collapsed;
            }
        }

        public object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
