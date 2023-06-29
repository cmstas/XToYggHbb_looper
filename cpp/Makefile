CC=g++
CFLAGS=-I.
OUT=main.exe 
OBJS=main.cc 

run_Hgg: $(OUT)

$(OUT): $(OBJS) ScanChain_Hgg.C
	$(CC) ${OBJS} -o $(OUT) `root-config --cflags --glibs` ../NanoCORE/NANO_CORE.so
	
clean:
	rm main.exe
