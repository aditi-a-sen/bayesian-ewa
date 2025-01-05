import pandas as pd

# Load the data into a Pandas DataFrame
data = pd.read_csv("patent_race_with_drug.csv")  # Replace "your_data.csv" with the actual file path

# Get unique participants
unique_participants = data["SID"].unique()

# Create an empty list to store the selected data
selected_data = []

# Iterate over each participant
for participant in unique_participants:
    # Filter data for the current participant
    participant_data = data[data["SID"] == participant]

    # Iterate over each session
    for session in participant_data["Session"].unique():
        # Filter data for the current session
        session_data = participant_data[participant_data["Session"] == session]

        # Iterate over each drug
        for drug in session_data["Drug"].unique():
            # Filter data for the current drug
            drug_data = session_data[session_data["Drug"] == drug]

            # Extract the first 5 and last 5 rounds
            first_5_rounds = drug_data.head(5)
            last_5_rounds = drug_data.tail(5)

            # Combine the first and last 5 rounds
            selected_data.append(pd.concat([first_5_rounds, last_5_rounds]))

# Combine the selected data into a single DataFrame
final_data = pd.concat(selected_data)

# Print or save the final data
print(final_data)
print("done!")
final_data.to_csv("selected_data.csv", index=False)