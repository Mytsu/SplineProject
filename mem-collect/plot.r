#
#Script gerado automaticamente por mem-collector.sh
#
dados <- read.table(file='coleta-memoria-11.dat');
png('coleta-memoria-11.png', width=1200);
plot(dados$V1, dados$V2, type="l", main="Uso de Memoria",
xlab="Tempo (1 pto/0.25 seg)", ylab="Mem (kb)");
points(dados$V1, dados$V2, col="red", pch=20);
dev.off();
