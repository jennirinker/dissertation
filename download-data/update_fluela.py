
if (__name__ == '__main__'):
    """ Do when called as script from command line.    
    """
    from URLDataDownload import updateDirectory

    # 20-Hz base URL and directory
    #   *** NO TRAILING SLASHES ***
    baseURL = 'http://wind.nrel.gov/MetData/135mData/' + \
                'FluelaTwr/20Hz/mat'
#    basedir = 'G:\\data\\fluela-high_freq'
    basedir = 'C:\\Users\\jrinker\\Documents\\temp_data\\fluela-high_freq'

    # update the directory, return list of err'd files
    errList = updateDirectory(baseURL,basedir)

    # print a comment if any errors
    if errList:
        print('\n{} file(s) '.format(len(errList)) \
            + 'with errors saved in errList.')

    # and that's it!
    print '\nScript complete.\n'

    
