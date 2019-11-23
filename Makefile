CC = g++
IDIR = include

CFLAGS = -I $(IDIR)

SRC = $(wildcard *.cpp)
NAME = a.out


$(NAME): $(SRC)
	@$(CC) -o $(NAME) $(CFLAGS) $(SRC)

.PHONY: clean
clean:
	@rm -f $(NAME)
