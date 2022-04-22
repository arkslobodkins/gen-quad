A1 = rand(40, 30);
[QA1,R] = qr(A1);

fileID = fopen('A1.txt','w');
fprintf(fileID, '%i %i ', size(A1, 1), size(A1, 2));
fprintf(fileID, '%.16e ', A1);
fclose(fileID);


fileID = fopen('Q1.txt','w');
fprintf(fileID, '%i %i ', size(QA1, 1), size(QA1, 2));
fprintf(fileID, '%.16e ', QA1);
fclose(fileID);