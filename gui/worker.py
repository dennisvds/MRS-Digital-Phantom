from PyQt5.QtCore import QThread, pyqtSignal

class Worker(QThread):
    finished = pyqtSignal(object)     # To send the final result
    message = pyqtSignal(str)          # To send text/log updates
    progress = pyqtSignal(int)         # To send progress (0-100%)

    def __init__(self, func, *args, **kwargs):
        super().__init__()
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def run(self):
        try:
            # Inject log_message and log_progress functions if needed
            if 'log_message' in self.func.__code__.co_varnames:
                self.kwargs['log_message'] = self.emit_message
            if 'log_progress' in self.func.__code__.co_varnames:
                self.kwargs['log_progress'] = self.emit_progress

            result = self.func(*self.args, **self.kwargs)
            self.finished.emit(result)
        except Exception as e:
            self.finished.emit(e)

    def emit_message(self, text):
        self.message.emit(text)

    def emit_progress(self, percent):
        self.progress.emit(percent)
