CC = mpic++
IDIR = include

CFLAGS = -I $(IDIR) -O2

SRC = $(wildcard *.cpp)
NAME = a.out


$(NAME): $(SRC)
	@$(CC) -o $(NAME) $(CFLAGS) $(SRC) -Wno-narrowing

.PHONY: clean
clean:
	@rm -f $(NAME) data/* runfile.* movie_runfile.* *.mp4 times.csv
