# Structured logs in `modkit motif search`

The debug logs can me emitted using the `--log-filepath` option to `modkit motif search`.
These logs are JSON-formatted lines that can be queried with a tool such as [jq](https://jqlang.org/) or [jaq](https://github.com/01mf02/jaq).
The top-level schema is:

| name        | description                                                       |
|-------------|-------------------------------------------------------------------|
| timestamp   | system time of the log message                                    |
| level       | log level, ERROR, WARN, INFO, or DEBUG                            |
| fields      | JSON object with more information about the event (details belos) |
| target      | Module logging the message (not usually very useful)              |
| filename    | Filename where the log message originated                         |
| line_number | Line in the file where the log message originated                 |

For example, to print out the INFO logs to the terminal:

```bash
cat ${log} | jq 'select(.level == "INFO") | .fields.message'
```

The `fields` object contains more information that can be useful for drilling down into what happened during a search. 

| name       | required | description                                                                                             |
|------------|----------|---------------------------------------------------------------------------------------------------------|
| message    | true     | The human-readable log message                                                                          |
| mod_code   | false    | The modification code (e.g. `a`) being worked on when this event was emitted                            |
| stage      | false    | The stage of the search algorithm, one of {`seeded`, `seedless`, `search`}                              |
| motif      | false    | The motif under consideration, this is the hunam-readable name e.g. `G[a]TC`                            |
| from_motif | false    | In the refinement step, this is the "input" motif, usually a seed, the same notation as `motif` is used |
| action     | false    | One of {`found`, `refined`, `discard`}                                                                  |
| require    | false    | Required value (only present when `action = discard`)                                                   |
| value      | false    | The value this motif has (only present when `action = discard`)                                         |

### Actions
- `found` A motif passes all criteria, and is kept.
  (It may later be decided that it was _re_-found and `discard`ed).
- `refined` When a motif "seed" is transformed into a motif that passes criteria, usually the step right before `found`, but connects this motif to the seed or `from_motif`
- `discard` When a motif is proposed by `refined`, but fails for some reason, such as it was already discovered

## Examples:

- Messages for a motif:
```bash
cat ${log} | jq 'select(.fields.motif == "TCG[a]") | .fields.message'
```
- Information for when a motif was refined:

```bash
cat ${log} | jq 'select(.fields.from_motif == "R[a]GC") | .fields'
```
- Check all of the seeds that were refined:

```bash
cat ${log} | jq 'select(.fields.stage == "seeded" and .fields.action == "found" and .fields.mod_code == "a") | .fields'
```

