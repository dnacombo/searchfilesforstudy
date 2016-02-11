function eegplugin_std_search( fig, try_strings, catch_strings)
plotmenu = findobj(fig, 'tag', 'EEGLAB');
studysub = findobj(plotmenu, 'label', 'Create study');
uimenu(studysub, 'label', 'Search files...', 'callback', ...
    [try_strings.no_check '[STUDY ALLEEG ]= pop_std_search( ALLEEG,EEG,CURRENTSET );eeglab redraw;' catch_strings.add_to_hist  ]);
