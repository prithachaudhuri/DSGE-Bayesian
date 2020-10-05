function [f]=sssolver(g, sigmas, sigmau, alphaa, lambda, xiH, xiL, gammas, gammau, eta, musw, muuw, mup, etan, etaw)
f1= g(1)-g(37)*mup-(1-alphaa)*(g(13)-g(40)*mup);
f2= g(2)-(1-alphaa)*g(14)-1;
f3= g(3)-(1-alphaa)*g(15);
f4= g(4)-(1-alphaa)*g(16);
f5= g(25)-g(43)*mup-log(1-alphaa)+alphaa*(g(13)-g(40)*mup);
f6= g(26)-1+alphaa*g(14);
f7= g(27)+alphaa*g(15);
f8= g(28)+alphaa*g(16);
f9= g(29)-g(44)*mup-musw-xiH-sigmas*(g(5)-g(38)*mup)-gammas*(g(17)-g(41)*mup);
f10= g(30)-sigmas*g(6)-gammas*g(18);
f11= g(31)-sigmas*g(7)-gammas*g(19)-1;
f12= g(32)-sigmas*g(8)-gammas*g(20);
f13= g(33)-g(45)*mup-muuw-sigmau*(g(9)-g(39)*mup)-gammau*(g(21)-g(42)*mup);
f14= g(34)-sigmau*g(10)-gammau*g(22);
f15= g(35)-sigmau*g(11)-gammau*g(23);
f16= g(36)-sigmau*g(12)-gammau*g(24)-1;
f17= g(9)-g(39)*mup-g(33)-g(21)+g(45)*mup+g(42)*mup;
f18= g(10)-g(34)-g(22);
f19= g(11)-g(35)-g(23);
f20= g(12)-g(36)-g(24);
f21= g(17)-g(41)*mup-g(21)+g(42)*mup+eta*(g(29)-g(44)*mup-g(33)+g(45)*mup);
f22= g(18)-g(22)+eta*(g(30)-g(34));
f23= g(19)-g(23)+eta*(g(31)-g(35));
f24= g(20)-g(24)+eta*(g(32)-g(36));
f25= g(2)-lambda*g(6)*exp(g(5)-g(1))-(1-lambda)*g(10)*exp(g(9)-g(1));
f26= g(3)-lambda*g(7)*exp(g(5)-g(1))-(1-lambda)*g(11)*exp(g(9)-g(1));
f27= g(4)-lambda*g(8)*exp(g(5)-g(1))-(1-lambda)*g(12)*exp(g(9)-g(1));
f28= lambda*exp(g(5))+(1-lambda)*exp(g(9))-exp(g(1));
f29= g(14)-lambda*g(18)*(exp(g(17)-g(13))^etan)-(1-lambda)*g(22)*(exp(g(21)-g(13))^etan);
f30= g(15)-lambda*g(19)*(exp(g(17)-g(13))^etan)-(1-lambda)*g(23)*(exp(g(21)-g(13))^etan);
f31= g(16)-lambda*g(20)*(exp(g(17)-g(13))^etan)-(1-lambda)*g(24)*(exp(g(21)-g(13))^etan);
f32= lambda*(exp(g(17))^etan)+(1-lambda)*(exp(g(21))^etan)-exp(g(13))^etan;
f33= g(26)-lambda*g(30)*(exp(g(29)-g(25))^etaw)-(1-lambda)*g(34)*(exp(g(33)-g(25))^etaw);
f34= g(27)-lambda*g(31)*(exp(g(29)-g(25))^etaw)-(1-lambda)*g(35)*(exp(g(33)-g(25))^etaw);
f35= g(28)-lambda*g(32)*(exp(g(29)-g(25))^etaw)-(1-lambda)*g(36)*(exp(g(33)-g(25))^etaw);
f36= lambda*(exp(g(29))^etaw)+(1-lambda)*(exp(g(33))^etaw)-exp(g(25))^etaw;
f37= g(37)-(1-alphaa)*g(40);
f38= g(43)+alphaa*g(40)+1;
f39= g(44)-sigmas*g(38)-gammas*g(41);
f40= g(45)-sigmau*g(39)-gammau*g(42);
f41= g(39)-g(45)-g(42);
f42= g(41)-g(42)+eta*(g(44)-g(45));
f43= g(37)-lambda*g(38)*exp(g(5)-g(1))-(1-lambda)*g(39)*exp(g(9)-g(1));
f44= g(40)-lambda*g(41)*(exp(g(17)-g(13))^etan)-(1-lambda)*g(42)*(exp(g(21)-g(13))^etan);
f45= g(43)-lambda*g(44)*(exp(g(29)-g(25))^etaw)-(1-lambda)*g(45)*(exp(g(33)-g(25))^etaw);

f=[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13;f14;f15;f16;f17;f18;f19;f20;f21;f22;f23;f24;f25;f26;f27;f28;f29;f30;f31;f32;f33;f34;f35;f36;f37;f38;f39;f40;f41;f42;f43;f44;f45];
end