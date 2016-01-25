function result =accelerometerError (co , acc)
 scaleMatrix = [ co(1)  co(4)  co(5);
                co(4)  co(2)  co(6);
                co(5)  co(6)  co(3);];
%  scaleMatrix = [ co(1)  0  0;
%                 0  co(2)  0;
%                 0  0 co(3);];
acc(:,1) = acc(:,1) - co(7);
acc(:,2) = acc(:,2) - co(8);
acc(:,3) = acc(:,3) - co(9);

result = (scaleMatrix*acc(:,1:3)')';

result = sqrt(sum((result.^2)')');

% 
% Error(:) = mean(Error)*length(Error);


