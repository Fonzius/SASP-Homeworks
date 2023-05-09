function plotFeatures3d(labels_unwrap, phi_unwrap)

m1 = labels_unwrap;
m1(m1 == 2 | m1 == 3) = 0;
k1 = m1.* phi_unwrap;

m2 = labels_unwrap;
m2(m2 == 1 | m2 == 3) = 0;
m2(m2 == 2) = 1;
k2 = m2.* phi_unwrap;

m3 = labels_unwrap;
m3(m3 == 1 | m3 == 2) = 0;
m3(m3 == 3) = 1;
k3 = m3.* phi_unwrap;

figure()
plot3(k1(:,1), k1(:,2), k1(:,3), '.', 'Color','r');
hold on
plot3(k2(:,1), k2(:,2), k2(:,3), '.', 'Color','g');
hold on
plot3(k3(:,1), k3(:,2), k3(:,3), '.', 'Color','b');
xlabel('A1');
ylabel('A2');
zlabel('P');
legend('f1', 'f2', 'f3');

end