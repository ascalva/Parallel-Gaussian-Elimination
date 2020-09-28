import matplotlib.pyplot as plt
import pandas            as pd
import seaborn           as sns

sns.set()

def plot_time_vs_N(df):

    # Create plot object.
    fig, axes = plt.subplots(dpi=70, figsize=(12, 8))

    # Setup plot.
    axes.set_title("Execution Time vs. Number of Elements (N) for Different Values of Processors",
                   fontsize = 14)

    axes.set_xlim([1, 4000])
    axes.set_ylim([0, 70])

    axes.set_xlabel("Number of elements (N)", fontsize = 12)
    axes.set_ylabel("Total Execution Time (s)", fontsize = 12)

    for attr in df.columns[1:] :
        axes.plot(df["N"], df[attr], marker='o', label = f"P = {int(attr[1:])}", alpha = 0.8)

    plt.legend()
    plt.savefig("time_v_N.png", dpi=100)
    plt.show()


def main():

    # Read in data into data frame.
    df        = pd.read_csv("./data.csv")

    plot_time_vs_N(df)

main()
