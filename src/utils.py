import os


def checkDir(directory):
    """ Analagous to mkdir -p directory from the command line
    """
    if not os.path.isdir(directory):
        print("Dir %s doesn't exist. Creating it" % (directory))
        try:
            os.makedirs(directory)
        except OSError:
            # if multiple parallel processes are running, they could be trying to create the same directory at the same time
            print("Dir %s already exists" % (directory))
