function a = prop_2state(x, p)
R = x(1);
P = x(2);

a = [p.kR;
     p.gR*R;
     p.kP*R;
     p.gP*P];
end