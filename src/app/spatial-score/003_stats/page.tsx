"use client";

import { useState } from "react";

const MIN_GROUPS = 2;
const MAX_GROUPS = 4;

type Group = { name: string; files: File[] };

const emptyGroup = (): Group => ({ name: "", files: [] });

export default function Boxplot() {
    const [groups, setGroups] = useState<Group[]>([emptyGroup(), emptyGroup()]);

    const addGroup = () => setGroups(g => [...g, emptyGroup()]);
    const removeGroup = (i: number) => setGroups(g => g.filter((_, idx) => idx !== i));

    const updateName = (i: number, name: string) =>
        setGroups(g => g.map((grp, idx) => idx === i ? { ...grp, name } : grp));

    const updateFiles = (i: number, files: FileList | null) =>
        setGroups(g => g.map((grp, idx) => idx === i ? { ...grp, files: files ? Array.from(files) : [] } : grp));

    return (
        <div>
            {groups.map((grp, i) => (
                <div key={i}>
                    <input
                        type="text"
                        value={grp.name}
                        onChange={e => updateName(i, e.target.value)}
                        placeholder="Group name"
                    />
                    <input
                        type="file"
                        multiple
                        onChange={e => updateFiles(i, e.target.files)}
                    />
                    {groups.length > MIN_GROUPS && (
                        <button onClick={() => removeGroup(i)}>-</button>
                    )}
                </div>
            ))}
            {groups.length < MAX_GROUPS && (
                <button onClick={addGroup}>+</button>
            )}
        </div>
    );
}
