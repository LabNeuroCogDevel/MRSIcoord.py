function scout3 = ShowCircle_test(scoutf,res,pl,pr,radius)
   scout3 = scoutf;
   % show the circle for 3
   for i=1:res
      for j=1:res
         % calculate the distance
         dist = round(sqrt((pl(1)-i)^2 + (pl(2)-j)^2));
         if (dist==round(radius/2.0))

            scout3(i,j,1)=1;
            scout3(i,j,2)=1;
            scout3(i,j,3)=0;

         end
      end
   end
   imshow(flipud(scout3),'InitialMagnification','fit');
   % AdjustImage_test(1,scout3,scout3);
   pause(0.5);
   for i=1:res
      for j=1:res
         dist = round(sqrt((pr(1)-i)^2 + (pr(2)-j)^2));
         if (dist==round(radius/2.0))

            scout3(i,j,1)=1;
            scout3(i,j,2)=1;
            scout3(i,j,3)=0;

         end        
      end


   end


   imshow(flipud(scout3),'InitialMagnification','fit');
   % AdjustImage_test(1,scout3,scout3);
   % pause(0.5);
