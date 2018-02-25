using System;
using System.Windows.Input;

namespace Processing
{
    public class DelegateCommand : ICommand
    {
        public Action<object> ExecuteAction { set; get; }
        public Func<object, bool> CanExecuteFunc { get; set; }

        public event EventHandler CanExecuteChanged;

        public DelegateCommand(Action<object> execute)
        {
            ExecuteAction = execute;
        }

        public bool CanExecute(object parameter)
        {
            if (CanExecuteFunc != null)
            {
                return CanExecuteFunc(parameter);
            }
            else
            {
                return true;
            }
        }

        public void Execute(object parameter)
        {
            ExecuteAction?.Invoke(parameter);
        }

        public void RaiseCanExecuteChanged()
        {
            CanExecuteChanged?.Invoke(this, EventArgs.Empty);
        }
    }
}
