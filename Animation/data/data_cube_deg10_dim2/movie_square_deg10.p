# Create movie with mencoder
ENCODER = system('which mencoder');
if (strlen(ENCODER)==0) print '=== mencoder not found ==='; exit
CMD = 'mencoder mf://@file-list -mf fps=1:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o square_deg10.avi'
system(CMD)
