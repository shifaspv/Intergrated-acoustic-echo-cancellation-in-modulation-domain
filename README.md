# Intergrated-acoustic-echo-cancellation-in-modulation-domain

In this repository included the MATLAB implementation of Modulation domain joint echo and noise supression technique proposed in the paper Cited below.

The echo signal was created on run with a specific impulse response selected by the user internally in the main script "AES_echonoise_modulation.m". In addition to echo, the near-end backround noise is also added into the near-end speech, which is Babble noise in the current script. User can create a .MAT file for their noise of wish and feed it into the model.

The main function is "AES_echonoise_modulation.m", in which the entire task of echo and noise estimation and cancellation has been done. The code has commented internally, reader can go and see the details inside.
All other functions are secondary, which are being called internally to the above function for different purposes (creating echo signal, framming speech, etc)


# Cite as
Jayakumar, E. P., PV Muhammed Shifas, and P. S. Sathidevi. "Integrated acoustic echo and noise suppression in modulation domain." International Journal of Speech Technology 19.3 (2016): 611-621.
