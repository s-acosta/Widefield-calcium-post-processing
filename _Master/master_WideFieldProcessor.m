% WideFieldProcessor Master

session_filename = 'WTR040_201205.tif';

%%

WF = WideFieldProcessor(session_filename);
WF.readStack()
WF.computeProperties()
WF.registerDFF()

    


