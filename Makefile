CC = mpic++
IDIR = include

CFLAGS = -I $(IDIR) -O2

SRC = $(wildcard *.cpp)
NAME = FDTD_sim


$(NAME): $(SRC)
	@$(CC) -o $(NAME) $(CFLAGS) $(SRC) -Wno-narrowing

.PHONY: clean
clean:
	@rm -f $(NAME) data/* 
