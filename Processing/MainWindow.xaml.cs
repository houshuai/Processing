using System.Windows;

namespace Processing
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            DataContext = new ViewModel();
            InitializeComponent();
        }

        private void EdgeDetectionMenu_Click(object sender, RoutedEventArgs e)
        {
            EdgeDetectionPopup.IsOpen = true;
        }

        private void SeparationMenu_Click(object sender, RoutedEventArgs e)
        {
            SeparationPopup.IsOpen = true;
        }

        private void EdgeDetectionStartBtn_Click(object sender, RoutedEventArgs e)
        {
            EdgeDetectionPopup.IsOpen = false;
        }

        private void SeparationStartBtn_Click(object sender, RoutedEventArgs e)
        {
            SeparationPopup.IsOpen = false;
        }
    }
}
