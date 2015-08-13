function runclassexamples()
% Run all CPLEX Class API examples

% ---------------------------------------------------------------------------
% File: runclassexamples.m
% Version 12.6.1
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
% Copyright IBM Corporation 2008, 2014. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
% ---------------------------------------------------------------------------

disp ('++lpex1,r');
lpex1 ('r');

disp ('++lpex1,c');
lpex1 ('c');

disp ('++lpex1,a');
lpex1 ('a');

disp ('++lpex2,afiro.mps');
lpex2 ('../../data/afiro.mps');

disp ('++lpex3');
lpex3();

disp ('++lpex6');
lpex6();

disp ('++lpex7,afiro.mps');
lpex7 ('../../data/afiro.mps');


% MIP
disp ('++mipex1,r');
mipex1 ('r');

disp ('++mipex1,c');
mipex1 ('c');

disp ('++mipex1,a');
mipex1 ('a');

disp ('++mipex2,mexample.mps');
mipex2 ('../../data/mexample.mps');

disp ('++mipex3,r');
mipex3 ('r');

disp ('++mipex3,c');
mipex3 ('c');

disp ('++mipex3,a');
mipex3 ('a');

disp ('++mipex4,noswot.mps');
mipex4 ('../../data/noswot.mps');


% QP and QCP
disp ('++qcpex1');
qcpex1;

disp ('++qcpdual');
qcpdual;

disp ('++qpex1');
qpex1;

disp ('++qpex2,qpex.lp');
qpex2 ('../../data/qpex.lp');

disp ('++indefqpex1');
indefqpex1;


% Others
disp ('++ miqpex1')
miqpex1;

disp ('++ cutstock')
cutstock;

disp ('++ rates');
rates;

disp ('++ facility');
facility;

disp ('++ blend');
blend;

disp ('++ fixcost1');
fixcost1;

disp ('++ diet');
diet ('r');

disp ('++ diet');
diet ('c');

disp ('++ diet');
diet ('i');

disp ('++ inout1');
inout1;

disp ('++ inout3');
inout3;

disp ('++ populate, ../../data/location.lp');
populate ('../../data/location.lp');

disp ('++ steel');
steel;

disp ('++ foodmanu');
foodmanu;

disp ('++ fixcost1');
fixcost1;

end

