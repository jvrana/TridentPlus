from pydent.browser import Browser
import os
import dill


class SavedBrowser(object):

    class RECORDMODES(object):

        ALWAYS = 'ALWAYS'         # record all browser changes upon exit
        ONCE = 'ONCE'  # record upon exit only if the browser file did not exist

    def __init__(self, filepath, session, record_mode=RECORDMODES.ONCE):
        self.filepath = filepath
        self.session = session
        self.browser = None
        self.record_mode = record_mode
        self.is_new_browser = False

    def load_browser(self):
        if os.path.isfile(self.filepath):
            with open(self.filepath, 'rb') as f:
                print('')
                browser = dill.load(f)
                self.is_new_browser = False
        else:
            browser = Browser(self.session)
            self.is_new_browser = True
        self.browser = browser
        return self.browser

    def save_browser(self):
        with open(self.filepath, 'wb') as f:
            dill.dump(self.browser, f)

    def __enter__(self):
        self.load_browser()
        return self.browser

    def __exit__(self, type, value, traceback):
        if self.record_mode == self.RECORDMODES.ALWAYS:
            self.save_browser()
        elif self.record_mode == self.RECORDMODES.ONCE and self.is_new_browser:
            self.save_browser()
