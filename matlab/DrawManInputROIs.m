function DrawManInputROIs(scout2)

   % Victor Yushmanov 04/08/2020
   % Quick and dirty, needs improvements
   % Derived from CalcPositionsThalamus.m in the original SVR1H2019 package
   % To be used when SVR1H2019 GUI is open and images loaded
   % Calls ShowCircle_test.m

   eh = findall(0,'Tag','RESscout');
   res = str2num(get(eh,'String'))

   % get the radius
   eh = findall(0,'Tag','Radius');
   radius = str2num(get(eh,'String'));

   eh = findall(0,'Tag','L1Row');
   posl(1,1) = str2num(get(eh,'String'));
   eh = findall(0,'Tag','L1Col');
   posl(1,2) = str2num(get(eh,'String'));
   eh = findall(0,'Tag','R1Row');
   posr(1,1) = str2num(get(eh,'String'));
   eh = findall(0,'Tag','R1Col');
   posr(1,2) = str2num(get(eh,'String'));

   %     note res +2 refers such that center point of image does not switch
   %     row entry flipped by res+2 bec MATLAB rows matrix starts on top (left, top
   %     is origin)
   %     col entry flipped by res+2 bec MATLAB cols matrix starts on left
   disp('these coords are in SID orientation but w transformed/cropped scout');
   posl(1,:) = res+2-posl(1,:);
   posr(1,:) = res+2-posr(1,:);

   ShowCircle_test(scout2,res,posl(1,:),posr(1,:),radius);

