CC = mpic++
IDIR = include

CFLAGS = -I $(IDIR) -O3

SRC = $(wildcard *.cpp)
NAME = a.out


$(NAME): $(SRC)
	@$(CC) -o $(NAME) $(CFLAGS) $(SRC)

.PHONY: clean
clean:
	@rm -f $(NAME) data/*
