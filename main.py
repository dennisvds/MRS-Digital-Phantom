import sys
from PyQt5.QtWidgets import QApplication
from gui.main_gui import MainWindow

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

def preprocess_metab_df():
    from utils.preprocess_df import create_metab_df
    create_metab_df(labels=['Background', 'WM', 'GM', 'CSF'], groups=['Healthy', 'Control'], 
                                                fraction_boundary=0.6, age_range=[18, 60], save=True)

if __name__ == '__main__':
    # Uncomment the following line to preprocess the metabolite DataFrame first
    # preprocess_metab_df()

    # Start the main application
    main()
    