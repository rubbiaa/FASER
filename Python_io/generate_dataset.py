import argparse
import ROOT
import numpy as np
import glob
import tqdm
from ROOT import TFile

def get_channel_xyz_from_id(ID, geom_detector):
    """
    Calculate the x, y, z coordinates from the hit ID and geometry detector.
    Args:
        ID (int): The hit ID.
        geom_detector: The geometry detector object.
    Returns:
        tuple: The x, y, z coordinates.
    """
    hittype = ID // 100000000000

    if hittype == 0:  # hit in scintillator
        ix = ID % 1000
        iy = (ID // 1000) % 1000
        iz = (ID // 1000000) % 1000
        ilayer = ID // 1000000000

        x = ix * geom_detector.fScintillatorVoxelSize - geom_detector.fScintillatorSizeX / 2.0
        y = iy * geom_detector.fScintillatorVoxelSize - geom_detector.fScintillatorSizeY / 2.0
        z = (ilayer * geom_detector.fSandwichLength + iz * geom_detector.fScintillatorVoxelSize
            - (geom_detector.NRep * geom_detector.fSandwichLength) / 2.0
            + geom_detector.fScintillatorVoxelSize * 2)
        return x, y, z

    elif hittype == 1:  # hit in Si tracker
        ix = ID % 10000
        iy = (ID // 10000) % 10000
        ilayer = (ID // 100000000) % 1000

        x = ix * geom_detector.fSiTrackerPixelSize - geom_detector.fScintillatorSizeX / 2.0
        y = iy * geom_detector.fSiTrackerPixelSize - geom_detector.fScintillatorSizeY / 2.0
        z = (ilayer * geom_detector.fSandwichLength + geom_detector.fSandwichLength
            - (geom_detector.NRep * geom_detector.fSandwichLength) / 2.0)
        return x, y, z

    else:
        print(f"TcalEvent::getChannelXYZfromID - hit of unknown type {hittype}")
        return 0, 0, 0

def main():
    """
    Main function to process ROOT files and convert them to Numpy format.
    """
    parser = argparse.ArgumentParser(description='Process ROOT files and convert to Numpy format.')
    parser.add_argument('--lib', required=True, help='Path to the shared object library.')
    parser.add_argument('--input', required=True, help='Directory containing ROOT files.')
    parser.add_argument('--output', required=True, help='Directory to save Numpy converted files.')
    
    args = parser.parse_args()

    # Load the shared object library
    ROOT.gSystem.Load(args.lib)
    files = glob.glob(args.input + "/*")
    print(f"Number of files: {len(files)}")

    # Loop through ROOT files
    t = tqdm.tqdm(enumerate(files), total=len(files), disable=False)
    for filenumber, filepath in t:
        f = TFile(filepath, "read")
        tree = f.Get("calEvent")
        
        # Check if the tree exists and can be accessed
        if not tree:
            print(f"File {filepath}: unable to access 'calEvent' tree.")
            continue

        tree.GetEntry(0)
        filename = filepath.split("/")[-1][18:-5]

        try:
            tracks = tree.tracks
            geom = tree.geom
        except AttributeError:
            print(f"File {filepath}: wrong format")
            continue

        all_hits_list = []

        # Process each track in the tree
        for track in tracks:
            track_id = track.ftrackID
            parent_id = track.fparentID
            primary_id = track.fprimaryID
            pdg = track.fPDG
            hits = track.fhitIDs
            energy_dep = track.fEnergyDeposits

            hit_info = np.zeros((len(hits), 8))

            for i, hit_id in enumerate(hits):
                x, y, z = get_channel_xyz_from_id(hit_id, geom)
                hit_info[i] = [track_id, parent_id, primary_id, pdg, x, y, z, energy_dep[i]]

            all_hits_list.append(hit_info)

        if all_hits_list:
            all_hits = np.concatenate(all_hits_list)
            np.savez_compressed(f"{args.output}/{filenumber}", hits=all_hits, filename=filename)

        del all_hits_list
        f.Close()

if __name__ == '__main__':
    main()

