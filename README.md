# srp-bigint

## How to run
1. Install [Docker](https://docs.docker.com/desktop/) and [Docker Compose](https://docs.docker.com/desktop/).
2. Download the source code.

```git```

3. Build the image, you only have to do this once.

```docker compose build```

4. Run the image. Results will be saved in the results directory.

```docker compose up```

5. To choose what to benchmark, edit the bottom of `testBigInt.cpp`, then re-run the image.