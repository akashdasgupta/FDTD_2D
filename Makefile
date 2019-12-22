CC = mpic++
IDIR = include

CFLAGS = -I $(IDIR) -O1

SRC = $(wildcard *.cpp)
NAME = a.out


$(NAME): $(SRC)
	@$(CC) -o $(NAME) $(CFLAGS) $(SRC)

.PHONY: clean
clean:
	@rm -f $(NAME) data/* runfile.* 
