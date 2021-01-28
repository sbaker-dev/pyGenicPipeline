import eel

FORM_VALUES = {}

@eel.expose
def set_genetic_path(path_text):
    """
    This will set the path to FORM_VALUES from a given path provided via path_text, if it cannot be used then return an
    error

    :param path_text: A loadable .psd path
    :type path_text: str

    :return: Nothing, just update the path or the error message based on the type
    """
    print(path_text)
    if len(path_text) == 0:
        eel.showStatus("You did not give a path!", "genetic_path_set", True)
    else:
        eel.showStatus("Set path!", "genetic_path_set", False)
        FORM_VALUES["set_genetic_path"] = path_text
        print(FORM_VALUES)


if __name__ == '__main__':
    eel.init("Web")

    eel.start("index.html")