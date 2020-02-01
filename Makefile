CC = mpic++
IDIR = include

CFLAGS = -I $(IDIR) -fopenmp -O3

SRC = $(wildcard *.cpp)
NAME = a.out


$(NAME): $(SRC)
	@$(CC) -o $(NAME) $(CFLAGS) $(SRC) -Wno-narrowing

.PHONY: clean
clean:
	@rm -f $(NAME) data/* runfile.* movie_runfile.* *.mp4
